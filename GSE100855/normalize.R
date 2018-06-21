library(DESeq2)
library(ggplot2)
library(vsn)  # meanSdPlot
library(sva)
library(mmalign)

out.fname <- filename("GSE100855");


pheno <- qread("samples.tsv");

in.pattern <- "\\.txt$";
in.fpaths <- list.files("raw", pattern=in.pattern, full.names=TRUE);
in.fnames <- list.files("raw", pattern=in.pattern, full.names=FALSE);
xs <- lapply(in.fpaths, read.table, header=FALSE, sep="\t", colClasses=c("character", "numeric"));

# verify that all gene names are the same
genes <- xs[[1]][,1];
stopifnot(unlist(lapply(xs, function(x) all(x[,1] == genes))))

# collect all counts together
x <- as.data.frame(lapply(xs, function(x) x[,2]));
samples <- gsub("(GSM[^_]+)_.*", "\\1", in.fnames);
colnames(x) <- samples;
rownames(x) <- genes;
x <- as.matrix(x);


#

dds <- DESeqDataSetFromMatrix(countData = x, colData = pheno, design = ~ treatment + cell);

# Fit DESeq model
dds <- DESeq(dds);

# list the coefficients
resultsNames(dds)

# partial knockout of protein expression
# (homozygous inframe deletion; hemizygous missense)
ko142 <- results(dds, name="cell_HeLa.KO142_vs_HeLa");
ko142[genes == "XPA", ]
summary(ko142)

# complete knockout of protein expression
# (homozygous frameshift deletion)
ko38 <- results(dds, name="cell_HeLa.KO38_vs_HeLa");
ko38[genes == "XPA", ]
summary(ko38)

# complemented with XPA cDNA
xp12ro.xpa <- results(dds, contrast=c("cell", "XP12RO.XPA", "XP12RO"));
xp12ro.xpa[genes == "XPA", ]
summary(xp12ro.xpa)

# complemented with XPA cDNA
xp20s.xpa <- results(dds, contrast=c("cell", "XP20S.XPA", "XP20S"));
xp20s.xpa[genes == "XPA", ]
summary(xp20s.xpa)

#

# shrunk coefficients are more useful for visualization
res <- lfcShrink(dds, coef="cell_HeLa.KO142_vs_HeLa", type="apeglm");
plotMA(res)
#idx <- identify(res$baseMean, res$log2FoldChange);
#genes[idx]

mcols(res)

#

d <- plotCounts(dds, gene=which(genes == "XPA"), intgroup="cell", returnData=TRUE);

qdraw(
	ggplot(d, aes(x=cell, y=count)) +
		geom_jitter(width=0.1) + scale_y_log10(breaks=c(50, 100, 500, 1000, 5000)) +
		theme_bw() + coord_flip()
	,
	height = 3,
	file = insert(out.fname, "expr-xpa", ext="pdf")
)

#d <- plotCounts(dds, gene=which(genes == "XPC"), intgroup="cell", returnData=TRUE);
d <- plotCounts(dds, gene=which(genes == "DDB2"), intgroup="cell", returnData=TRUE);
d <- plotCounts(dds, gene=which(genes == "ERCC4"), intgroup="cell", returnData=TRUE);
ggplot(d, aes(x=cell, y=count)) +
	geom_jitter(width=0.1) + scale_y_log10(breaks=c(50, 100, 500, 1000, 5000)) +
	theme_bw() + coord_flip()

#

# last 12 samples have ~2x the total counts cmopared to the first 36
cbind(pheno, colSums(x))

lib.sizes <- colSums(x);
qdraw(
	ggplot(data.frame(pheno, lib_size = colSums(x)),
		aes(x=factor(description, levels=rev(description)), y=lib_size)) +
		geom_bar(stat="identity") + theme_bw() + coord_flip() +
		ylab("Library size") + xlab("")
	,
	height = 6,
	file = insert(out.fname, "lib-size", ext="pdf")
)

keep <- rowSums(counts(dds)) >= 10;
dds.f <- dds[keep, ];

rld <- rlog(dds, blind=TRUE);
y <- assay(rld);

hist(log(x + 1), breaks=1000)
hist(y, breaks=1000)
hist(y[y > 0], breaks=1000)
meanSdPlot(y)

rld.f <- rlog(dds.f, blind=TRUE);
y.f <- assay(rld.f);

summary(as.numeric(y.f))
hist(y.f, breaks=1000)
hist(y.f[y.f > 0], breaks=1000)
hist(y.f, breaks=100)
meanSdPlot(y.f)

vsd <- vst(dds.f, blind=TRUE);
y.v <- assay(vsd);

summary(as.numeric(y.v))
hist(y.v, breaks=100)
hist(y.v[y.v > min(y.v)], breaks=1000)
meanSdPlot(y.v)

hist((y.f - min(y.f)) - (y.v - min(y.v)))

ntd <- normTransform(dds.f);
y.n <- assay(ntd);

summary(as.numeric(y.n))
hist(y.n[y.n > min(y.n)], breaks=1000)
hist(y.n[y.n > min(y.n)], breaks=100)
meanSdPlot(y.n)

hist(log(rowSums(counts(dds)) + 1), breaks=1000)


boxplot(y)
y2 <- y;
y2[y <= 0] <- 0;
boxplot(y2)
boxplot(y.f)
boxplot(y.v)
boxplot(y.n)

cor(as.numeric(y.f), as.numeric(y.v))
cor(as.numeric(y.f), as.numeric(y.n))

# VST has a squashing effect on low values
qdraw(
	{
		smoothScatter(as.numeric(y.f), as.numeric(y.v))
		abline(a=0, b=1, col="orangered")
	},
	file = insert(out.fname, c("sm-scatter", "rlg-vst"), ext="pdf")
)

# Negative values for y.f should probably be thresholded
qdraw(
	{
		smoothScatter(as.numeric(y.f), as.numeric(y.n), xlab="Regularized log transform", ylab="Log transform")
		abline(a=0, b=1, col="orangered")
	},
	file = insert(out.fname, c("sm-scatter", "rlg-log"), ext="pdf")
)

# VST has a squashing effect on low values
qdraw(
	{
		smoothScatter(as.numeric(y.n), as.numeric(y.v), xlab="Log transform", ylab="Variance-stabilizing transform")
		abline(a=0, b=1, col="orangered")
	},
	file = insert(out.fname, c("sm-scatter", "log-vst"), ext="pdf")
)

#

y2 <- y.v;

# cell is the variable of interest
# use no adjustment variable
mod <- model.matrix(~cell, data=pheno);
mod0 <- model.matrix(~1, data=pheno);

#n.sv <- num.sv(y2, mod, method="be");

sv <- sva(y2, mod, mod0, n.sv=1);

barplot(sv$sv[,1])
#barplot(sv$sv[,2])
#barplot(sv$sv[,3])

fsv <- fsva(y2, mod, sv, y2);

summary(as.numeric(y2 - fsv$db))
hist(as.numeric(y2 - fsv$db), breaks=100)

y.sva <- fsv$db;

hist(y.sva, breaks=1000)


#

#library(pheatmap)

#idx <- order(apply(y, 1, sd), decreasing=TRUE)[1:1000];
#pheatmap(y[idx, ], show_rownames=FALSE, annotation_col=pheno)

#y.pca <- pca(y);
#y.pca <- pca(y.f);
y.pca <- pca(y.v);
#y.pca <- pca(y.n);

qdraw(
	# PC1 is dominated by parent cell lineage
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment)) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("pca", "pc-1-2"), ext="pdf")
)

qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=2:3) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("pca", "pc-2-3"), ext="pdf")
)

pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=3:4) +
	scale_shape_manual(values=c(1, 19))

# PC4 is dominated by batch? (last 12 samples vs. first 36)
qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=4:5) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("pca", "pc-4-5"), ext="pdf")
)

# PC6 may be XPA response
# PC7 may be retinoic acid response
qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=6:7) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("pca", "pc-6-7"), ext="pdf")
)

# PC8 and PC9 may also bee retinoic acid response
qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=8:9) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("pca", "pc-8-9"), ext="pdf")
)


#

#d <- dist(1 - cor(y.f));
d <- dist(1 - cor(y.v));
dm <- as.matrix(d);
rownames(dm) <- colnames(dm) <- pheno$description;

qdraw(
	{
		pheatmap(dm, clustering_distance_rows=d, clustering_distance_cols=d, cex=0.8)
	},
	width = 8, height = 8,
	file = insert(out.fname, tag="cor-heatmap", ext="pdf")
);

#

# y.sva has been adjusted for the hidden batch (first 36 vs. last 12)
y.pca <- pca(y.sva);

qdraw(
	# PC1 is dominated by parent cell lineage
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment)) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("sva", "pca", "pc-1-2"), ext="pdf")
)

qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=2:3) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("sva", "pca", "pc-2-3"), ext="pdf")
)

pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=3:4) +
	scale_shape_manual(values=c(1, 19))

# PC4 is dominated by batch? (last 12 samples vs. first 36)
# (not after SVA)
qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=4:5) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("sva", "pca", "pc-4-5"), ext="pdf")
)

# PC 6 may be XPA effect
qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=6:7) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("sva", "pca", "pc-6-7"), ext="pdf")
)

# look for most negative values, since samples with negative PC6 values
# are XPA+
pc <- y.pca$U["PC6", ];
sort(pc, decreasing=FALSE)[1:100]
hist(pc, breaks=100);

qqnorm(pc); qqline(pc);

qdraw(
	pca_plot(y.pca, pheno=pheno, mapping=aes(colour=cell, size=xpa, shape=treatment), dims=8:9) +
		scale_shape_manual(values=c(1, 19))
	,
	width = 6.5,
	file = insert(out.fname, c("sva", "pca", "pc-8-9"), ext="pdf")
)


#

qwrite(dds, insert(out.fname, "deseq-data-set", ext="rds"));

qwrite(x, insert(out.fname, "expr-counts", ext="mtx"));
qwrite(y.f, insert(out.fname, "expr-rlg", ext="mtx"));
qwrite(y.v, insert(out.fname, "expr-vst", ext="mtx"));
qwrite(y.n, insert(out.fname, "expr-log", ext="mtx"));
qwrite(y.sva, insert(out.fname, "expr-vst-sva", ext="mtx"));

