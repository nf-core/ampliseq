workflow PARSE_INPUT {
    take:
    input // folder
    single_end
    multiple_sequencing_runs
    extension

    main:
    // Folder input

    //Check folders in folder when multiple_sequencing_runs
    folders = multiple_sequencing_runs ? "/*" : ""
    error_message = "\nCannot find any reads matching: \"${input}${folders}${extension}\"\n"
    error_message += "Please revise the input folder (\"--input_folder\"): \"${input}\"\n"
    error_message += "and the input file pattern (\"--extension\"): \"${extension}\"\n"
    error_message += "*Please note: Path needs to be enclosed in quotes!*\n"
    error_message += multiple_sequencing_runs ? "If you do not have multiple sequencing runs, please do not use \"--multiple_sequencing_runs\"!\n" : "If you have multiple sequencing runs, please add \"--multiple_sequencing_runs\"!\n"
    error_message += "In any case, please consult the pipeline documentation.\n"
    if ( single_end ) {
        //Get files - single end
        Channel
            .fromPath( input + folders + extension )
            .ifEmpty { error("${error_message}") }
            .map { read ->
                    def meta = [:]
                    meta.sample           = read.baseName.toString().indexOf("_") != -1 ? read.baseName.toString().take(read.baseName.toString().indexOf("_")) : read.baseName
                    meta.single_end   = single_end.toBoolean()
                    meta.run          = multiple_sequencing_runs ? read.take(read.findLastIndexOf{"/"})[-1] : "1"
                    [ meta, read ] }
            .set { ch_reads }
    } else {
        //Get files - paired end
        Channel
            .fromFilePairs( input + folders + extension, size: 2 )
            .ifEmpty { error("${error_message}") }
            .map { name, reads ->
                    def meta = [:]
                    meta.sample           = name.toString().indexOf("_") != -1 ? name.toString().take(name.toString().indexOf("_")) : name
                    meta.single_end   = single_end.toBoolean()
                    meta.run          = multiple_sequencing_runs ? reads[0].take(reads[0].findLastIndexOf{"/"})[-1] : "1"
                    [ meta, reads ] }
            .set { ch_reads }
    }
    if (multiple_sequencing_runs) {
        //Get folder information
        ch_reads
            .flatMap { meta, reads -> [ meta.run ] }
            .unique()
            .set { ch_folders }
        //Report folders with sequencing files
        ch_folders
            .collect()
            .subscribe {
                String folders = it.toString().replace("[", "").replace("]","")
                log.info "\nFound the folder(s) \"$folders\" containing sequencing read files matching \"${extension}\" in \"${input}\".\n" }
        //Stop if folder count is 1 and multiple_sequencing_runs
        ch_folders
            .count()
            .subscribe { if ( it == 1 ) error("Found only one folder with read data but \"--multiple_sequencing_runs\" was specified. Please review data input.") }
    }

    //Check whether all sampleID = meta.sample are unique
    ch_reads
        .map { meta, reads -> [ meta.sample ] }
        .toList()
        .subscribe {
            if( it.size() != it.unique().size() ) {
                ids = it.take(10);
                error("Please review data input, sample IDs are not unique! First IDs are $ids")
            }
        }

    //Check that no dots "." are in sampleID
    ch_reads
        .map { meta, reads -> meta.sample }
        .subscribe { if ( "$it".contains(".") ) error("Please review data input, sampleIDs may not contain dots, but \"$it\" does.") }

    //Check that sampleIDs do not start with a number when using metadata (sampleID gets X prepended by R and metadata wont match any more!)
    ch_reads
        .map { meta, reads -> meta.sample }
        .subscribe { if ( params.metadata && "$it"[0].isNumber() ) error("Please review data input, sampleIDs may not start with a number, but \"$it\" does. The pipeline unintentionally modifies such strings and the metadata will not match any more.") }

    emit:
    reads   = ch_reads
}
