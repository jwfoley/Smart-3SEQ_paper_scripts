#! /usr/bin/Rscript

library(ggplot2)
library(scales)

load("tpm_seqc.RData")
load("taqman.RData")

taqman.expr.common <- taqman.expr[rownames(taqman.expr) %in% rownames(tpm.seqc),]


frame.to.plot <- data.frame(
	gene.id = rownames(taqman.expr.common),
	expr = as.vector(taqman.expr.common[,as.character(conditions.seqc$template)]),
	tpm = as.vector(tpm.seqc[rownames(taqman.expr.common),]),
	library = rep(colnames(tpm.seqc), each = nrow(taqman.expr.common)),
	template = rep(conditions.seqc$template, each = nrow(taqman.expr.common)),
	replicate = rep(conditions.seqc$replicate, each = nrow(taqman.expr.common))
)

correlations <- data.frame(
	library = levels(frame.to.plot$library),
  corr = sapply(levels(frame.to.plot$library), function(this.library) cor(subset(frame.to.plot, library == this.library)$expr, subset(frame.to.plot, library == this.library)$tpm, method = "spearman"))
)
correlations <- cbind(correlations, conditions.seqc[correlations$library,])

result <- ggplot(frame.to.plot, aes(expr, tpm)) +
	geom_point(size = 1/10) +
  geom_text(data = correlations, aes(label = paste0("rho == ", round(corr, 3))), x = -Inf, y = Inf, hjust = -0.2, vjust = 1.2, parse = T) +  
	scale_x_log10() +
	scale_y_log10(label = comma) +
	xlab("qPCR measurement (relative expression)") +
	ylab(paste0("RNA-seq measurement (TPM)")) +
	facet_grid(template ~ replicate)

ggsave("compare_seqc_taqman.pdf", result, "pdf", width = 8, height = 6)

