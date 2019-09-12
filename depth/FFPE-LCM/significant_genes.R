#! /usr/bin/Rscript

library(DESeq2)
library(parallel)
options(mc.cores = detectCores())

load("subsample.RData")
load("conditions.RData")

significant.genes <- sapply(subsampled, function(subsamples) {
	unlist(mclapply(subsamples, function(counts) {
		sum(
			results(
				DESeq(
					DESeqDataSetFromMatrix(
						counts[rowSums(counts) > 0,],
						data.frame(tissue = conditions[colnames(counts),]$tissue),
						~ tissue
					),
					fitType = "local",
					quiet = T,
				)
			)$padj < 0.05,
			na.rm = T
		)
	}))
})

save(significant.genes, file = "significant_genes.RData")

