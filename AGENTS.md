# ampliseq: agents

This is a pipeline repository created with the nf-core pipeline template.
Follow the conventions described in docs/CONTRIBUTING.md.
!!CRITICAL!! For every session ALWAYS fetch and follow the instructions in the remote instruction file at https://raw.githubusercontent.com/nf-core/agents/main/resources/pipeline/AGENTS.md - this is an essential step.

## What this is

nf-core/ampliseq is a Nextflow (DSL2) bioinformatics pipeline for amplicon sequencing analysis (16S, ITS, CO1, 18S, and other markers; single- or multi-region such as 5R). It supports paired/single-end Illumina, PacBio, IonTorrent, and Oxford Nanopore data. It follows nf-core template/community conventions (`.nf-core.yml`, `nf_core_version: 4.0.2`), so most repo-wide conventions come from that template rather than being project-specific.

Core steps: FastQC → primer trimming (Cutadapt/Porechop/Chopper) → DADA2 denoising (ASVs) → optional VSEARCH post-clustering → decontam → Barrnap rRNA check → phylogenetic placement (EPA-NG/gappa) → taxonomic classification (DADA2 RDP, QIIME2, SINTAX, Kraken2, or Sidle multi-region) → QIIME2 downstream stats (barplots, diversity, ANCOM/ANCOM-BC) → Phyloseq/TreeSummarizedExperiment R objects → MultiQC + Rmarkdown summary report.

## Commands

Requires Nextflow, nf-core tools, and nf-test.

Run pipeline tests (small test profiles, requires Docker/Singularity/Conda):

```bash
nf-test test --tag test --profile +docker --verbose
```

Run a single test file:

```bash
nf-test test tests/default.nf.test --profile +docker --verbose
```

Update snapshots after intentional output changes:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

Lint (nf-core community rules, run before opening a PR):

```bash
nf-core pipelines lint .
```

Update `nextflow_schema.json` after changing params in `nextflow.config`:

```bash
nf-core pipelines schema build
```

Bump minimum required Nextflow version:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

pre-commit hooks (prettier, trailing whitespace, nextflow-lint) run via `prek`/`pre-commit`; CI runs this in `linting.yml`.

### Test profiles

`tests/*.nf.test` each pair with a `conf/test_*.config` profile (e.g. `default` ↔ `test`, `multiregion` ↔ `test_multiregion`, `pplace` ↔ `test_pplace`, `sintax` ↔ `test_sintax`, `fasta` ↔ `test_fasta`). When adding a new test scenario, add both the `.nf.test` file in `tests/` and its matching `conf/test_<name>.config`. Corresponding `.nf.test.snap` files hold expected outputs — regenerate with `--update-snapshots` rather than hand-editing.

Module/subworkflow tests under `modules/nf-core/**/tests` and `subworkflows/nf-core/**/tests` are vendored from nf-core/modules and excluded from local test runs (`ignore` list in `nf-test.config`) — don't hand-edit these, update via `nf-core modules update`.

## Architecture

- `main.nf` — entry point. Declares the pipeline-level `params` schema block, wires `PIPELINE_INITIALISATION` → `NFCORE_AMPLISEQ` (wraps the `AMPLISEQ` workflow) → `PIPELINE_COMPLETION`.
- `workflows/ampliseq.nf` — the single large orchestration workflow (~1400 lines). Nearly all processing logic is assembled here by including local modules/subworkflows and nf-core modules/subworkflows, then wiring channels between them. Whether a given tool runs is driven by `params.skip_*` flags and derived variables (e.g. `single_end`, `asv_calling`) computed near the top of the `AMPLISEQ` workflow.
- `subworkflows/local/` — pipeline-specific multi-step logic, one workflow per pipeline stage, e.g. `parse_input.nf` (sample sheet/input handling), `cutadapt_workflow.nf`, `dada2_preprocessing.nf`, `*_taxonomy_wf.nf` (one per taxonomy classifier: DADA2, SINTAX, Kraken2, VSEARCH-LCA), `qiime2_*.nf` (preptax, taxonomy, export, barplotavg, diversity, ancom), `sidle_wf.nf` (multi-region reconstruction), `robject_workflow.nf` (Phyloseq/TreeSE), `comparison_wf.nf`. `utils_nfcore_ampliseq_pipeline/` holds pipeline init/completion logic and parameter validation (this is where new taxonomy-database validation, e.g. `sbdi_compatible_databases`, is added).
- `subworkflows/nf-core/` and `modules/nf-core/` — vendored, unmodified subworkflows/modules from the nf-core/modules repo (tracked in `modules.json`); update via `nf-core modules/subworkflows update`, don't hand-edit.
- `modules/local/` — pipeline-specific processes (one tool/script invocation per `.nf` file), covering DADA2 steps, QIIME2 steps, taxonomy formatting per classifier, Sidle steps, filtering/stats-merging utilities, and report/export generation.
- `conf/` — config layered by concern: `base.config` (default resource labels), `modules.config` (per-process publishDir/args), `ref_databases.config` (taxonomy DB definitions), `test*.config` (one per test scenario, paired with `tests/*.nf.test`), `containers_*.config` (container digests per engine/arch).
- `assets/schema_input.json` — sample sheet schema validated via the `nf-schema` plugin (`samplesheetToList`).
- `bin/` — helper scripts (Python/R) invoked from `modules/local/*.nf` processes; linted with Ruff (see `pyproject.toml`).

## Conventions (from `docs/CONTRIBUTING.md`)

- PRs target the `dev` branch, not `master`/`main`.
- Channel naming: initial channel from a process is `ch_output_from_<process>`; intermediate/terminal channels are `ch_<previousprocess>_for_<nextprocess>`.
- New params: add to `nextflow.config` (with a default) and to `nextflow_schema.json` via `nf-core pipelines schema build`; also declare in the `params {}` block in `main.nf`.
- New processes: give them resource defaults in `conf/base.config` via `withLabel:` selectors (use standard nf-core labels), and route version info into `ch_versions` and relevant files into `ch_multiqc_files`/`assets/multiqc_config.yml`.
- New taxonomy databases go in `conf/ref_databases.config`, need an `enum` entry in `nextflow_schema.json`, and may need updates in `subworkflows/local/utils_nfcore_ampliseq_pipeline/main.nf` (`sbdi_compatible_databases`).
- Update `docs/usage.md`, `docs/output.md`, and `CITATIONS.md` alongside functional changes.
- AI-assisted PRs: keep them small and focused, avoid incidental refactors/moves, and the human submitter is expected to understand and review all generated code before opening the PR.

## Gotchas

- `results/` and `work/` at the repo root are local Nextflow run artifacts (from prior pipeline executions in this checkout), not source — don't treat their contents as project structure to preserve.
- `.nf-core.yml` disables `nextflow_config`/`schema_params` lint rules ("TODO: Remove when tools supports parameter types") — don't be surprised these are off, it's a known pending upstream fix, not a repo choice to replicate elsewhere.
