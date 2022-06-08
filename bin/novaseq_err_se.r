#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
    stop("Usage: novaseq_err_se.r <model> <run_id>")
}

model <- args[1]
run_id <- args[2]

suppressPackageStartupMessages(library(dada2))

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

errF = readRDS(model)

#monotone decreasing
sink(file = paste0(run_id, ".md.err.log"))
make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))

errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
errF.md.full <- errF
errF.md.full$err_out <- errF.md
saveRDS(errF.md.full, paste0(run_id, ".md.err.rds"))
sink(file = NULL)

pdf(paste0(run_id, ".md.err.pdf"))
plotErrors(errF.md.full, nominalQ = TRUE)
dev.off()

sink(file = paste0(run_id, ".md.err.convergence.txt"))
dada2:::checkConvergence(errF.md.full)
sink(file = NULL)