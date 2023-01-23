#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
    stop("Usage: novaseq_err_pe.r <fw_model> <rv_model> <run_id>")
}

fw_model <- args[1]
rv_model <- args[2]
run_id <- args[3]

suppressPackageStartupMessages(library(dada2))

errF = readRDS(fw_model)
errR = readRDS(rv_model)

#monotone decreasing
make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))

errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
colnames(errF.md) <- colnames(errF$err_out)
errF.md.full <- errF
errF.md.full$err_out <- errF.md
saveRDS(errF.md.full, paste0(run_id, "_1.md.err.rds"))

errR.md <- t(apply(getErrors(errR), 1, make.monotone.decreasing))
colnames(errR.md) <- colnames(errR$err_out)
errR.md.full <- errR
errR.md.full$err_out <- errR.md
saveRDS(errR.md.full, paste0(run_id, "_2.md.err.rds"))

pdf(paste0(run_id, "_1.md.err.pdf"))
plotErrors(errF.md.full, nominalQ = TRUE)
dev.off()
svg(paste0(run_id, "_1.md.err.svg"))
plotErrors(errF.md.full, nominalQ = TRUE)
dev.off()

pdf(paste0(run_id, "_2.md.err.pdf"))
plotErrors(errR.md.full, nominalQ = TRUE)
dev.off()
svg(paste0(run_id, "_2.md.err.svg"))
plotErrors(errR.md.full, nominalQ = TRUE)
dev.off()

sink(file = paste0(run_id, "_1.md.err.convergence.txt"))
dada2:::checkConvergence(errF.md.full)
sink(file = NULL)

sink(file = paste0(run_id, "_2.md.err.convergence.txt"))
dada2:::checkConvergence(errR.md.full)
sink(file = NULL)
