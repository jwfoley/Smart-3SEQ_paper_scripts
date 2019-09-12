#! /usr/bin/Rscript

library(parallel)
options(mc.cores = detectCores())

increment <- 5E5

subsample <- function(counts, p) sapply(counts, function(count) rbinom(1, count, p))


load("counts.RData")
total.reads <- read.table("total_reads", col.names = c("library", "reads"))
total.reads.lookup <- total.reads$reads
names(total.reads.lookup) <- total.reads$library

subsampled <- lapply(levels(conditions$amount), function(this.amount) {
	these.libraries <- rownames(subset(conditions, amount == this.amount))
	counts.filtered <- read.counts[,these.libraries]
	targets <- increment * 1:floor(min(total.reads.lookup[these.libraries]) / increment)
	result <- mclapply(targets, function(target) {
		sapply(these.libraries, function(this.library) subsample(counts.filtered[,this.library], target / total.reads.lookup[this.library]))
	})
	names(result) <- sapply(targets, format, scientific = F)
	result
})
names(subsampled) <- levels(conditions$amount)

save(subsampled, conditions, increment, file = "subsample.RData")

