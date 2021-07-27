// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def parse_samplesheet(LinkedHashMap row, single_end) {
    //Check if manifest contains column sampleID  & forwardReads
    if (row.sampleID == null || row.forwardReads == null) {
        exit 1, "ERROR: Please check input samplesheet -> Column 'sampleID' and 'forwardReads' are required but not detected."
    }
    //Check if manifest contains a column for reverse reads
    if (row.reverseReads == null && !single_end) {
        exit 1, "ERROR: Please check input samplesheet -> Column 'reverseReads' is missing. In case you do have only single ended reads, please specify '--single_end', '--pacbio', or '--iontorrent'."
    }
    //read meta info
    def meta = [:]
    meta.id           = row.sampleID
    meta.single_end   = single_end.toBoolean()
    meta.run          = row.run == null ? "1" : row.run
    //read data info
    def array = []
    if (!file(row.forwardReads).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Forward read FastQ file does not exist!\n${row.forwardReads}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.forwardReads) ] ]
    } else {
        if (!file(row.reverseReads).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Reverse read FastQ file does not exist!\n${row.reverseReads}"
        }
        array = [ meta, [ file(row.forwardReads), file(row.reverseReads) ] ]
    }
    return array
}