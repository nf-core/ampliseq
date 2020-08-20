#!/usr/bin/env Rscript
#########################################
# dada2_denoise_pacbio
#
# Run DADA2 denoising, using sequence files in a given folder
# and reading error model from given rds file.
#
# Jeanette Tångrot

suppressPackageStartupMessages(library(optparse))

VERSION = 1.0

# Get arguments
option_list = list(
  make_option(
    c('--filterDir'), type='character', default='dada2_filtered',
    help='Directory containing quality filtered reads to estimate sequence errors from, default "dada2_filtered".'
  ),
  make_option(
    c('--errModel'), type='character', default='err.rds',
    help='R RDS file with calculated error model, as generated by DADA2 learnErrors. Default: "err.rds".'
  ),
make_option(
    c('--prefix'), type='character', default='./',
    help='Prefix for names of generated files. Default: "./"'
  ),
make_option(
    c('--pool'), type='character', default='TRUE',
    help='Whether to pool together all samples prior to sample inference. Possible options: TRUE, pseudo, FALSE. Default: TRUE'
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
  write(sprintf("dada2_denoise_pacbio.r version %s, DADA2 version %s", VERSION, packageVersion('dada2')), stderr())
  q('no', 0)
}

# Check options
if ( ! file_test("-d", opt$filterDir) ) {
   stop( sprintf("Cannot find folder with filtered files: %s. See help (-h)\n",opt$filterDir) )
}

if ( ! file.exists(opt$errModel) ) {
   stop(sprintf("Cannot find %s. See help (-h).\n",opt$errModel))
}

if ( opt$pool == "TRUE" || opt$pool == "T") {
   opt$pool = TRUE
} else if
 ( opt$pool == "FALSE" || opt$pool == "F") {
   opt$pool = FALSE
} else if ( is.character(opt$pool) && opt$pool != "pseudo" ) {
   stop(sprintf("Invalid pool argument for dada2 denoising. See help (-h).\n",opt$errModel))
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

logmsg( sprintf( "Sequence denoising with DADA2." ) )

# Load DADA2 library here, to avoid --help and --version taking so long
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

# Write versions used
write(sprintf("%s\nDADA2: %s / Rcpp: %s / RcppParallel: %s", R.version.string, packageVersion('dada2'), packageVersion('Rcpp'), packageVersion('RcppParallel') ),file="")

# Read error model from file
err = readRDS(opt$errModel)

# Dereplicate identical reads
files=list.files(opt$filterDir,full.names=T)
derep <- derepFastq(files, verbose = opt$verbose)

# Denoising, save rds for dada2 object
dd <- dada(derep, err=err, multithread=T, pool=opt$pool)
saveRDS(dd,sprintf('%sdd.rds', opt$prefix))

logmsg(sprintf("Finished denoising"))
