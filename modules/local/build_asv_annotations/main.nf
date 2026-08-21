process BUILD_ASV_ANNOTATIONS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    path(barrnap_summary,    stageAs: 'barrnap/*')
    path(decontam_details,   stageAs: 'decontam/*')
    path(notcontam_details,  stageAs: 'notcontam/*')
    path(ssu_pre,            stageAs: 'ssu_pre/*')
    path(ssu_post,           stageAs: 'ssu_post/*')
    path(len_asv_pre,        stageAs: 'len_asv_pre/*')
    path(len_asv_post,       stageAs: 'len_asv_post/*')
    path(codons_pre,         stageAs: 'codons_pre/*')
    path(codons_post,        stageAs: 'codons_post/*')
    path(len_itsx_pre,       stageAs: 'len_itsx_pre/*')
    path(len_itsx_post,      stageAs: 'len_itsx_post/*')
    path(accept_pre,         stageAs: 'accept_pre/*')
    path(accept_post,        stageAs: 'accept_post/*')

    output:
    path("ASV_annotations.tsv"), emit: tsv
    path "versions.yml"         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Groovy-side presence checks -- each source is only wired in when the upstream step that
    // produces it actually ran (see workflows/ampliseq.nf); NULL here means "not run", not "no data"
    def barrnap_r      = barrnap_summary   ? "'${barrnap_summary}'"   : "NULL"
    def decontam_r     = decontam_details  ? "'${decontam_details}'" : "NULL"
    def notcontam_r    = notcontam_details ? "'${notcontam_details}'" : "NULL"
    def ssu_pre_r       = (ssu_pre       && ssu_post)      ? "'${ssu_pre}'"       : "NULL"
    def ssu_post_r      = (ssu_pre       && ssu_post)      ? "'${ssu_post}'"      : "NULL"
    def len_asv_pre_r   = (len_asv_pre   && len_asv_post)  ? "'${len_asv_pre}'"   : "NULL"
    def len_asv_post_r  = (len_asv_pre   && len_asv_post)  ? "'${len_asv_post}'"  : "NULL"
    def codons_pre_r    = (codons_pre    && codons_post)   ? "'${codons_pre}'"    : "NULL"
    def codons_post_r   = (codons_pre    && codons_post)   ? "'${codons_post}'"   : "NULL"
    def len_itsx_pre_r  = (len_itsx_pre  && len_itsx_post) ? "'${len_itsx_pre}'"  : "NULL"
    def len_itsx_post_r = (len_itsx_pre  && len_itsx_post) ? "'${len_itsx_post}'" : "NULL"
    def accept_pre_r    = (accept_pre    && accept_post)   ? "'${accept_pre}'"    : "NULL"
    def accept_post_r   = (accept_pre    && accept_post)   ? "'${accept_post}'"   : "NULL"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(readr)
    })

    barrnap_summary   <- $barrnap_r
    decontam_details  <- $decontam_r
    notcontam_details <- $notcontam_r
    ssu_pre           <- $ssu_pre_r
    ssu_post          <- $ssu_post_r
    len_asv_pre       <- $len_asv_pre_r
    len_asv_post      <- $len_asv_post_r
    codons_pre        <- $codons_pre_r
    codons_post       <- $codons_post_r
    len_itsx_pre      <- $len_itsx_pre_r
    len_itsx_post     <- $len_itsx_post_r
    accept_pre        <- $accept_pre_r
    accept_post       <- $accept_post_r

    # Barrnap: single winning-domain label per ASV, not the raw e-values (those are available in
    # barrnap's own native output). Argmin over non-NA e-value columns only -- NA (not some arbitrary
    # domain) when none are significant. NA here is deliberately ambiguous with "barrnap wasn't run"
    # (see docs/output.md); not disambiguated further, per confirmed design decision.
    chr_long <- tibble(asv_id = character(), field = character(), value = character())
    if (!is.null(barrnap_summary)) {
        evals <- read_tsv(barrnap_summary, show_col_types = FALSE)
        eval_cols <- setdiff(colnames(evals), c("ASV_ID", "eval_method"))
        label <- evals |>
            rowwise() |>
            mutate(barrnap_domain = {
                vals <- c_across(all_of(eval_cols))
                names(vals) <- sub("_eval\$", "", eval_cols)
                vals <- vals[!is.na(vals)]
                if (length(vals) == 0) NA_character_ else names(vals)[which.min(vals)]
            }) |>
            ungroup() |>
            transmute(asv_id = ASV_ID, field = "barrnap_domain", value = barrnap_domain)
        chr_long <- bind_rows(chr_long, label)
    }

    # Decontam and per-filter pass/fail: booleans only, kept in a separate (logical-typed) long
    # table from the character one above -- pivot_wider() has no type-guessing of its own, so
    # stacking different-typed sources into one generic value column would let bind_rows() silently
    # coerce booleans to 0/1 or strings, which then bakes the wrong column type into the Parquet
    # output. Splitting by type before stacking keeps every column's real type intact.
    lgl_long <- tibble(asv_id = character(), field = character(), value = logical())

    add_bool <- function(long, path, id_col, bool_col, field_name) {
        if (is.null(path)) return(long)
        df <- read_tsv(path, show_col_types = FALSE)
        bind_rows(long, tibble(asv_id = df[[id_col]], field = field_name, value = as.logical(df[[bool_col]])))
    }
    lgl_long <- lgl_long |>
        add_bool(decontam_details, "ID", "contaminant", "decontam_contaminant") |>
        add_bool(notcontam_details, "ID", "not.contaminant", "decontam_not_contaminant")

    # Filter pass/fail booleans: none of FILTER_SSU/FILTER_LEN/FILTER_CODONS emit a per-ASV boolean
    # today, only the surviving table -- derive it by diffing ASV_ID membership between the table as
    # it was immediately before this specific filter ran and immediately after (positional first
    # column, since the ID column's exact name isn't consistent across every filter's own output).
    add_filter_pass <- function(long, pre, post, field_name) {
        if (is.null(pre) || is.null(post)) return(long)
        pre_ids  <- read_tsv(pre,  show_col_types = FALSE)[[1]]
        post_ids <- read_tsv(post, show_col_types = FALSE)[[1]]
        bind_rows(long, tibble(asv_id = pre_ids, field = field_name, value = pre_ids %in% post_ids))
    }
    lgl_long <- lgl_long |>
        add_filter_pass(ssu_pre, ssu_post, "passed_ssu_filter") |>
        add_filter_pass(len_asv_pre, len_asv_post, "passed_length_filter_asv") |>
        add_filter_pass(codons_pre, codons_post, "passed_codon_filter") |>
        add_filter_pass(len_itsx_pre, len_itsx_post, "passed_length_filter_itsx") |>
        # Recapitulates the whole standard filtering chain (decontam through the ITSx-region length
        # filter) in one column, rather than combining the individual passed_*/decontam_* columns
        # above -- those need careful NA-vs-absent and decontam-polarity handling to combine
        # correctly (an ASV removed by an early filter is NA, not FALSE, in every later filter's own
        # column), which a direct start-to-end table diff sidesteps entirely.
        add_filter_pass(accept_pre, accept_post, "ampliseq_accept")

    wide_chr <- chr_long |> pivot_wider(names_from = field, values_from = value)
    wide_lgl <- lgl_long |> pivot_wider(names_from = field, values_from = value)
    if (nrow(wide_chr) == 0) wide_chr <- tibble(asv_id = character())
    if (nrow(wide_lgl) == 0) wide_lgl <- tibble(asv_id = character())

    annotations <- full_join(wide_chr, wide_lgl, by = "asv_id")
    write_tsv(annotations, "ASV_annotations.tsv")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major", "minor")], collapse = ".")),
            paste0("    r-dplyr: ", packageVersion("dplyr"))
        ),
        "versions.yml"
    )
    """

    stub:
    """
    echo -e "asv_id" > ASV_annotations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
