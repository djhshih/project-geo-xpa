library(affy)
library(frma)
library(io)
library(ggplot2)
library(reshape2)
library(mmalign)              # pca

library(hgu133plus2frmavecs)  # barcode
library(hgu133plus2.db)       # annotation


out.fname <- filename("GSE70818");


# Input ######################################################################

pheno <- qread("samples.tsv");
rownames(pheno) <- pheno$sample_id;

x <- ReadAffy(celfile.path="raw")

# populate phenoData
cel.names <- sampleNames(x);
sampleNames(x) <- sub("_.*\\.CEL", "", cel.names);
phenoData(x) <- AnnotatedDataFrame(
	pheno[match(sampleNames(x), pheno$sample_id), ]
);


# Normalize ##################################################################

# normalize data using fRMA
x.n <- frma(x);

# verify that all samples have median GNUSE < 1.25
# GNUSE = 1.25 indicates that precision is on average 25% worse than
# the typical array (of the large training data set)
gnuse.stats <- as.data.frame(t(GNUSE(x.n, type="stats")));
gnuse.stats <- data.frame(sample_id = rownames(gnuse.stats), gnuse.stats, check.names=FALSE);
qwrite(gnuse.stats, insert(out.fname, tag="gnuse-stats", ext="tsv"));

# plot GNUSE values
gnuse <- GNUSE(x.n, type="values");
qdraw(
	{
		par(mar=c(8, 6, 2, 2))
		boxplot(gnuse, outline=FALSE, las=2, ylab="GNUSE")
	},
	insert(out.fname, tag="gnuse", ext="pdf")
);

e <- exprs(x.n);

e.df <- melt(e, varnames=c("probe_id", "sample_id"));

qdraw(
	ggplot(e.df, aes(x=value)) +
		geom_density() + theme_bw() +
		facet_grid(sample_id ~ .) +
		xlab("fRMA normalized expression")
	,
	width = 4, height = 12,
	file = insert(out.fname, tag=c("density", "frma"), ext="pdf")
);

e.pca <- pca(e);

pca_plot(e.pca, pheno=pheno, mapping=aes(shape=cell, colour=factor(uv_treatment)));

# expression z-score using the barcode algorithm
# mu and sigma estimated from left-most mode across large dataset
# (unexpressed distribution)
# reduces the probe specific effect
z <- barcode(x.n, output="z-score");

z.df <- melt(z, varnames=c("probe_id", "sample_id"));

qdraw(
	ggplot(z.df, aes(x=value)) +
		geom_density() + theme_bw() +
		facet_grid(sample_id ~ .) +
		xlab("Expression barocode z-score")
	,
	width = 4, height = 12,
	file = insert(out.fname, tag=c("density", "barcode-z"), ext="pdf")
);

z.pca <- pca(z);

qdraw(
	pca_plot(z.pca, pheno=pheno, mapping=aes(shape=cell, colour=factor(time), size=factor(uv_treatment))) +
		scale_shape_manual(values=c(1, 16, 17))
	,
	width = 6.5,
	file = insert(out.fname, tag=c("pca", "barcode-z"), ext="pdf")
)


annot.db <- hgu133plus2.db;

probes <- rownames(z);
genes <- mapIds(annot.db, keys=probes, column="SYMBOL", keytype="PROBEID");
ensembls <- mapIds(annot.db, keys=probes, column="ENSEMBL", keytype="PROBEID");
entrezs <- mapIds(annot.db, keys=probes, column="ENTREZID", keytype="PROBEID");

# ensemble mapping rate is lower in the annotation database
mean(is.na(ensembls))
mean(is.na(genes))
mean(is.na(entrezs))


# Explore ####################################################################

# investigate how to aggregate probes to genes

m <- rowMeans(z);

sort(table(genes[!is.na(genes)]))
#probes.sel.idx <- which(genes == "YME1L1");
probes.sel.idx <- which(genes == "MALAT1");

plot(m[probes.sel.idx])

z.sel <- z[probes.sel.idx, ];
z.sel.df <- melt(z.sel, varnames=c("probe_id", "sample_id"));

ggplot(z.sel.df, aes(x=sample_id, y=value)) +
	geom_bar(stat="identity") +
	facet_grid(probe_id ~ .) + theme_bw()

z.sel.sd <- apply(z.sel, 1, sd);
z.sel.mean <- apply(z.sel, 1, mean);

z.sel.cmean <- colMeans(z.sel);
z.sel.cmean.df <- data.frame(sample_id=names(z.sel.cmean), value=z.sel.cmean)

ggplot(z.sel.cmean.df, aes(x=sample_id, y=value)) +
	geom_bar(stat="identity") + theme_bw()

# taking the probe with maximum signal => saturated probe, no variability
# taking the probe with minimum signal => weak binding probe
# taking the probe with highest sd => weak binding probe

# taking the mean is probably safest

map_probe_to_gene <- function(x, annot.db) {
	probes <- rownames(x);
	genes <- mapIds(annot.db, keys=probes, column="SYMBOL", keytype="PROBEID");
	x.g <- aggregate(x, list(genes), mean);
	rownames(x.g) <- x.g[,1];
	x.g <- x.g[,-1];
	as.matrix(x.g)
}


# Output #####################################################################

e.g <- map_probe_to_gene(e, annot.db);
z.g <- map_probe_to_gene(z, annot.db);

qwrite(e, insert(out.fname, tag=c("frma-expr", "probe"), ext="mtx"));
qwrite(e.g, insert(out.fname, tag=c("frma-expr", "gene"), ext="mtx"))

qwrite(z, insert(out.fname, tag=c("barcode-z", "probe"), ext="mtx"));
qwrite(z.g, insert(out.fname, tag=c("barcode-z", "gene"), ext="mtx"))

