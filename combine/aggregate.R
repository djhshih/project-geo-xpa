library(io)
library(sva)
library(dplyr)
library(mmalign)
library(vsn)
library(preprocessCore)

out.fname <- filename("geo-xpa");

expr.fnames <- c(
	#GSE100855 = "../GSE100855/GSE100855_expr-rlg.mtx",
	GSE100855 = "../GSE100855/GSE100855_expr-vst-sva.mtx",
	GSE55484 = "../GSE55484/GSE55484_bgrma-log-qn-expr_gene.mtx",
	GSE70818 = "../GSE70818/GSE70818_barcode-z_gene.mtx"
);

pheno.fnames <- c(
	GSE100855 = "../GSE100855/samples.tsv",
	GSE55484 = "../GSE55484/samples.tsv",
	GSE70818 = "../GSE70818/samples.tsv"
);

xs <- lapply(expr.fnames, qread); 
phenos <- mapply(
	function(fname, series) {
		x <- qread(fname);
		data.frame(study=series, x)
	},
	pheno.fnames, names(pheno.fnames)
);

gene.lists <- lapply(xs, rownames);


ner.genes <- c("XPA", "ERCC3", "XPC", "ERCC2", "DDB2", "ERCC4", "ERCC5", "POLH", "ERCC1", "ERCC6", "DDB1", "ERCC8", "UVSSA");

genes <- Reduce(intersect, gene.lists);

ner.genes %in% genes

x <- Reduce(
	function(x, y) {
		cbind(x[genes, ], y[genes, ])
	},
	xs
);

pheno <- Reduce(full_join, phenos);
pheno$uv_treatment[is.na(pheno$uv_treament)] <- 0;
pheno$time[is.na(pheno$uv_treament)] <- 0;

stopifnot(colnames(x) == pheno$sample_id)

#

boxplot(x)
meanSdPlot(x)

#x.v <- justvsn(x);

#boxplot(x.v)
#meanSdPlot(x.v)
#hist(x.v, breaks=1000)

x.q <- normalize.quantiles(x);
#x.q <- normalize.quantiles.robust(x, n.remove=0.3*ncol(x), use.log2=FALSE);
colnames(x.q) <- colnames(x);
rownames(x.q) <- rownames(x);

boxplot(x.q)
hist(x.q, breaks=1000)
meanSdPlot(x.q)

#x.qv <- justvsn(x.q);

#boxplot(x.qv)
#hist(x.qv, breaks=1000)
#meanSdPlot(x.qv)

#

mod.cb <- model.matrix(~1, data = pheno);
#x.cb <- ComBat(x.q, batch=pheno$study, mod=mod.cb, par.prior=TRUE, prior.plots=TRUE);
x.cb <- ComBat(x, batch=pheno$study, mod=mod.cb, par.prior=TRUE, prior.plots=TRUE);
#x.cb <- log(ComBat(exp(x), batch=pheno$study, mod=mod.cb, par.prior=TRUE,
#prior.plots=TRUE));
dev.off();

boxplot(x.cb)
hist(x.cb, breaks=200)
meanSdPlot(x.cb)

# VSN didn't help flatten sd vs. mean plot much
# need to exponentiate first because vsn does log-like transform
#x.cb.v <- justvsn(exp(x.cb));
#boxplot(x.cb.v)
#hist(x.cb.v, breaks=1000)
#meanSdPlot(x.cb.v)
#hist(x.cb.v - x.cb, breaks=100);


#

d <- data.frame(pheno, XPA=x.cb["XPA", ]);

qdraw(
	ggplot(d, aes(x=cell, y=XPA, colour=study)) +
		geom_jitter(width=0.1) +
		theme_bw() + coord_flip()
	,
	height = 3,
	file = insert(out.fname, "expr-xpa", ext="pdf")
)

#

x.pca <- pca(x);
#x.pca <- pca(x.q);

# first two PCs are dominated by study or platform specific effects
qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=study, colour=cell, size=xpa))
	,
	width = 7, height = 6,
	file = insert(out.fname, c("pca", "pc-1-2"), ext="pdf")
)

# PCs 3 and 4 show some XPA effect confounded with cell line effect
qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=study, colour=cell, size=xpa), dims=3:4)
	,
	width = 7, height = 6,
	file = insert(out.fname, c("pca", "pc-3-4"), ext="pdf")
)

# confounded effects
qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=study, colour=cell, size=xpa), dims=5:6)
	,
	width = 7, height = 6,
	file = insert(out.fname, c("pca", "pc-5-6"), ext="pdf")
)


x.pca <- pca(x.cb);

qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=study, colour=cell, size=xpa))
	,
	width = 7, height = 6,
	file = insert(out.fname, c("cb", "pca", "pc-1-2"), ext="pdf")
)

qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=study, colour=cell, size=xpa), dims=3:4)
	,
	width = 7, height = 6,
	file = insert(out.fname, c("cb", "pca", "pc-3-4"), ext="pdf")
)

qdraw(
	pca_plot(x.pca, pheno=pheno, mapping=aes(shape=study, colour=cell, size=xpa), dims=5:6)
	,
	width = 7, height = 6,
	file = insert(out.fname, c("cb", "pca", "pc-5-6"), ext="pdf")
)

#

qwrite(x.cb, insert(out.fname, "expr", ext="mtx"));
qwrite(x.cb, insert(out.fname, "expr", ext="rds"));
qwrite(pheno, insert(out.fname, "sample-info", ext="tsv"));

