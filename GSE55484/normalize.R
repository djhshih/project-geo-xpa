library(io)
library(lumi)
library(mmalign)
library(ggplot2)

library(illuminaHumanv4.db)


out.fname <- filename("GSE55484");

x.orig <- qread("raw/GSE55484_non-normalized.tsv");
pheno <- qread("samples.tsv");

# remove p-value detection columns
x <- x.orig[, grep("Pval", colnames(x.orig), invert=TRUE)];
# extract illumina id as rownames
rownames(x) <- x[,1];
x <- x[,-1];

data.frame(
	colnames(x),
	pheno
)

# NB  Due to lack of proper annotation, assume for now that
# samples in pheno data.frame and expression data are in the same order

colnames(x) <- pheno$sample_id;
x <- as.matrix(x);


# Normalize ##################################################################

# Since we will be merging this data with data from other platforms,
# we aim aim to reduce bias.

# Use RMA background correct to reduce bias

summary(as.numeric(x))
hist(log(x), breaks=100)

x.b <- lumiB(x, method="bgAdjust.affy");
#x.b <- lumiB(x, method="none");

summary(as.numeric(x.b))
summary(as.numeric(x.b - x))
hist(x.b - x, breaks=50)
hist(log(x.b), breaks=100)

boxplot(log(x.b), las=2, outline=FALSE)

x.t <- log2(x.b);

#x.n <- lumiN(x.t, method="rsn");
x.n <- lumiN(x.t, method="quantile");

summary(as.numeric(x.n))
summary(as.numeric(x.n - x.t))
hist(x.n, breaks=100)

boxplot(x.n, las=2, outline=FALSE)

# We cannot assess quality control by lumiQ because input data does not
# contain bead summary statistics (only expression and detection p-value
# provided).
# Hopefully, authors already discarded bad arrays...


# Diagnosis ##################################################################

x.pca <- pca(x.n);

# cell line effect is mixed together with xpa effect in the 
# first few principal components
# fourth component shows no biological signal of interest

qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=parental, colour=factor(xpa)))
	,
	width = 6,
	file = insert(out.fname, tag=c("pca", "qn-expr", "pc-1-2"), ext="pdf")
);

qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=parental, colour=factor(xpa)), dims=2:3)
	,
	width = 6,
	file = insert(out.fname, tag=c("pca", "qn-expr", "pc-2-3"), ext="pdf")
);

qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=parental, colour=factor(xpa)), dims=3:4)
	,
	width = 6,
	file = insert(out.fname, tag=c("pca", "qn-expr", "pc-3-4"), ext="pdf")
);

map_probe_to_gene <- function(x, annot.db) {
	probes <- rownames(x);
	genes <- mapIds(annot.db, keys=probes, column="SYMBOL", keytype="PROBEID");

	# aggregate probes to gene by mean
	x.g <- aggregate(x, list(genes), mean);
	rownames(x.g) <- x.g[,1];
	x.g <- x.g[,-1];

	as.matrix(x.g)
}

annot.db <- illuminaHumanv4.db;
x.n.g <- map_probe_to_gene(x.n, annot.db);

hist(x.n.g, breaks=100)


# Explore ####################################################################

# check biological positive control (XPA)

y <- data.frame(pheno, XPA = x.n.g["XPA", ]);

ggplot(y, aes(x=cell, y=XPA)) +
	geom_jitter(stat="identity", width=0.1) + theme_bw()

# however, XP20S-XPA does not seem to overexpress XPA by much

xpa.probes <- mapIds(annot.db, keys="XPA", keytype="SYMBOL", column="PROBEID");


y.r <- data.frame(pheno, XPA = x[xpa.probes[1], ]);

ggplot(y.r, aes(x=cell, y=XPA)) +
	geom_jitter(stat="identity", width=0.1) + theme_bw()

hist(log(x), breaks=100)
x.l.z <- scale(log(x));
hist(x.l.z, breaks=100)

y.z <- data.frame(pheno, XPA = x.l.z[xpa.probes[1], ]);

ggplot(y.z, aes(x=cell, y=XPA)) +
	geom_jitter(stat="identity", width=0.1) + theme_bw()


# Output #####################################################################

qwrite(x.n, insert(out.fname, tag=c("bgrma-log-qn-expr", "probe"), ext="mtx"));
qwrite(x.n.g, insert(out.fname, tag=c("bgrma-log-qn-expr", "gene"), ext="mtx"));

