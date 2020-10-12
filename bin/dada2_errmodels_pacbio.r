#!/usr/bin/env Rscript
#########################################
# dada2_errmodels_pacbio
#
# DADA2 sequence error estimation based on files in a given folder
# Assumes PacBio reads, i.e. uses errorEstimationFunction = PacBioErrfun
#
# Author: Jeanette Tångrot (jeanette.tangrot@nbis.se), Daniel Lundin

suppressPackageStartupMessages(library(optparse))

VERSION = 1.0

# Get arguments
option_list = list(
  make_option(
    c('--filterDir'), type='character', default='dada2_filtered',
    help='Directory containing quality filtered reads to estimate sequence errors from, default "dada2_filtered".'
  ),
  make_option(
    c('--prefix'), type='character', default='./',
    help='Prefix for name of rds file with DADA2 error model. Can include a path to another folder. Default: "./"'
  ),
  make_option(
    c('--nbases'), type='character', default=1e8,
    help='Minimum number of total bases to use for error estimation, please see DADA2 documentation for details. Default: 1e8'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE,
    help="Print progress messages."
  ),
  make_option(
    c("--version"), action="store_true", default=FALSE,
    help="Print version of script and DADA2 library."
  )
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$version ) {
  write(sprintf("dada2_errmodels_pacbio.r version %s, DADA2 version %s", VERSION, packageVersion('dada2')), stderr())
  q('no', 0)
}

# Check options
if ( ! file_test("-d", opt$filterDir) ) {
   stop( sprintf("Cannot find folder with filtered files: %s. See help (-h)\n",opt$filterDir) )
}

# Function for log messages
logmsg = function(msg, llevel='INFO') {
  if ( opt$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}

logmsg( sprintf( "Sequence error estimation with DADA2 learnErrors and PacBioErrFun." ) )

# Load DADA2 library here, to avoid --help and --version taking so long
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

# Do the error estimation, save rds for error profile
files=list.files(opt$filterDir,full.names=T)
logmsg( sprintf("Using files: %s", files))
err <- learnErrors(files, errorEstimationFunction=PacBioErrfun, multithread=TRUE, randomize=FALSE, verbose=opt$verbose, nbases=as.double(opt$nbases))
saveRDS(err,sprintf('%serr.rds', opt$prefix))

logmsg(sprintf("Finished error estimation"))
