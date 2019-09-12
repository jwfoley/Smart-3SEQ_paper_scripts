#! /usr/bin/Rscript

library(ggplot2)
library(scales)

total.reads <- read.table("ffpe-lcm_total_reads", as.is = T)
total.reads.lookup <- total.reads[,2]
names(total.reads.lookup) <- total.reads[,1]

star.logs <- read.table("ffpe-lcm_aggregate_logs.tsv", sep = "\t", header = T, as.is = T)

star.results <- data.frame(
	run = sub("\\.align\\.log$", "", star.logs$library),
	library = sub("\\..*", "", star.logs$library),
	uniquely.aligned = star.logs$Uniquely.mapped.reads.number,
	tissue = simplify2array(strsplit(star.logs$library, "_"))[1,],
	amount = simplify2array(strsplit(star.logs$library, "_"))[2,],
	length = as.integer(sub("^trunc", "", sapply(star.logs$library, function(x) strsplit(x, "\\.")[[1]][2])))
)

star.results$percent.uniquely.aligned <- star.results$uniquely.aligned / total.reads.lookup[as.character(star.results$library)]

result.plot <- ggplot(star.results, aes(length, percent.uniquely.aligned, group = library, color = amount)) +
	geom_line() +
	geom_vline(xintercept = c(51, 76), col = "darkgray", linetype = "dashed") +
	scale_y_continuous(label = percent) +
	xlab("read length (nt)") +
	ylab("uniquely aligned reads")

ggsave("read_length_ffpe-lcm.pdf", result.plot, width = 8, height = 4)

