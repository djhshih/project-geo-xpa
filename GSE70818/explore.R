library(io);
library(dplyr);
library(reshape2);
library(ggplot2);

z.g <- qread("GSE70818_barcode-z_gene.mtx");
pheno <- qread("samples.tsv");

ner.genes <- c("XPA", "ERCC3", "XPC", "ERCC2", "DDB2", "ERCC4", "ERCC5", "POLH", "ERCC1", "ERCC6", "DDB1", "ERCC8", "UVSSA");

z.ner <- melt(z.g[ner.genes, ], varnames=c("gene", "sample_id"), value.name="expression");

y <- left_join(pheno, z.ner);
y$uv_treatment <- factor(y$uv_treatment);

ggplot(y, aes(x=description, y=expression)) +
	geom_bar(stat="identity") + theme_bw() +
	facet_grid(gene ~ uv_treatment, scales="free_x") +
	theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplot(y, aes(x=description, y=expression)) +
	geom_bar(stat="identity") + theme_bw() +
	facet_grid(gene ~ time, scales="free_x") +
	theme(axis.text.x = element_text(angle = 45, hjust=1))

