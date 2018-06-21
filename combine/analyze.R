library(io)
library(limma)

x <- qread("geo-xpa_expr.rds");
pheno <- qread("geo-xpa_sample-info.tsv");

idx <- pheno$study == "GSE70818";
pheno$parental <- as.character(pheno$parental);
pheno$parental[idx] <- "NC";    # assume that all cells from that study are homogeneous
pheno$parental <- as.factor(pheno$parental);

pheno$xpa <- relevel(pheno$xpa, "+");

mod <- model.matrix(~ 0 + parental:xpa, pheno);
groups <- paste0(
	gsub("parental([A-Za-z0-9]+).*", "\\1", colnames(mod)),
	"_",
	ifelse(gsub(".*xpa(.+)", "\\1", colnames(mod)) == "+", "pos", "neg")
);
renamed <- data.frame(old=colnames(mod), new=groups);
colnames(mod) <- groups;

fit <- lmFit(x, mod);

n <- length(groups) / 2;
contrast.strs <- paste0(groups[(n+1):(2*n)], "-", groups[1:n]);
cmat <- makeContrasts(contrasts=contrast.strs, levels=groups);

fit2 <- contrasts.fit(fit, cmat);
fit2 <- eBayes(fit2);

groups[1]
topTable(fit2, coef=1)

groups[2]
topTable(fit2, coef=2)

groups[3]
topTable(fit2, coef=3)

groups[4]
topTable(fit2, coef=4)

groups[5]
topTable(fit2, coef=5)

ts <- fit2$t;
hist(ts, breaks=100)

t.m <- rowMeans(ts);
hist(t.m, breaks=100)

t.m.s <- sort(t.m);
t.m.s[1:100]
rev(t.m.s)[1:100]
