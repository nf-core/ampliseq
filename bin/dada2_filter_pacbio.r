#!/usr/bin/env Rscript
#########################################
# dada2_filter_pacbio
#
# Uses the filterAndTrim function in the dada2 package
# to filter long amplicon sequences generated with PacBio.
#
# Author: Jeanette Tångrot (jeanette.tangrot@nbis.se), Daniel Lundin

suppressPackageStartupMessages(library(optparse))

VERSION = 1.0

# Get arguments
option_list = list(
  make_option(
    c('--infile'), type='character', default='',
    help='Manifest file listing sample names and paths to files to filter and trim. No default.'
  ),
  make_option(
    c('--filterDir'), type='character', default='dada2_filtered',
    help='Directory for quality filtered reads, default "dada2_filtered". Will be created if it does not exist.'
  ),
  make_option(
    c('--stats'), type='character', default='filter_stats.tsv',
    help='File for writing filtering information. Default: "filter_stats.tsv".'
  ),
  make_option(
    c('--maxEE'), type='integer', default=-1,
    help='Maximum number of expected errors in sequence, default Inf.'
  ),
  make_option(
    c('--truncLen'), type='integer', default=0,
    help='Truncate sequence after truncLen bases, reads shorter than this are discarded. Default: 0 (no truncation).'
  ),
  make_option(
    c('--truncQ'), type='integer', default=2,
    help='truncQ option in filterAndTrim(). Default 2.'
  ),
  make_option(
    c('--minLen'), type='integer', default=20,
    help='Remove reads shorter than minLen, after trimming and truncation. Default: 20.'
  ),
  make_option(
    c('--maxLen'), type='integer', default=-1,
    help='Remove reads longer than maxLen, before trimming and truncation.Default: Inf.'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE,
    help="Print progress messages."
  ),
  make_option(
    c("--version"), action="store_true", default=FALSE,
    help="Print version of this script and of DADA2 library."
  )
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$version ) {
  write(sprintf("dada2_filter_pacbio.r version %s, DADA2 version %s", VERSION, packageVersion('dada2')), stderr())
  q('no', 0)
}

# Check options
if ( ! file.exists(opt$infile) ) {
stop(sprintf("Cannot find %s. See help (-h).\n",opt$infile))
}
if ( ! file_test("-d", opt$filterDir) ) { dir.create(opt$filterDir) }

if ( opt$maxEE < 0 ) { opt$maxEE = Inf }
if ( opt$maxLen < 0 ) { opt$maxLen = Inf }

# Function for log messages
logmsg = function(msg, llevel='INFO') {
  if ( opt$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

# Load DADA2 library here, to avoid --help and --version taking so long
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

logmsg( sprintf( "Read quality filtering with filterAndTrim. Options used:\n  maxN: %d, maxEE: %s, truncQ: %d, minLen: %d, maxLen: %s",
    0, opt$maxEE, opt$truncQ, opt$minLen, opt$maxLen )
)

# Do the filtering/trimming
input <- read.table(opt$infile, header = TRUE, sep = ",", colClasses = "character")

filt <- file.path(opt$filterDir, basename(input$absolute.filepath))
res_filt <- filterAndTrim(input$absolute.filepath, filt, maxN = 0, maxEE = opt$maxEE, truncQ = opt$truncQ, truncLen = opt$truncLen, minLen = opt$minLen, maxLen = opt$maxLen, compress = T, multithread = T, verbose = opt$verbose)

input["file"] <- basename(input$absolute.filepath)
output <- merge(input,res_filt, by.x="file", by.y ="row.names")

# Write filtering stats to file opt$stats
write.table( output, file = opt$stats, sep = "\t", row.names = FALSE, quote = FALSE)

logmsg(sprintf("Done filtering, output in %s", opt$filterDir))
