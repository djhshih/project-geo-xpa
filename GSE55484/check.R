library(GEOquery)
library(ggplot2)

x <- getGEO("GSE55484")[[1]];

# XPA is assayed by probe ILMN_1787591
# ILMN_1787591 has sequence TGAGGACAAGATACCAAGGCAAACCCTAGATTGGGGTAGAGGGAAAAGGG
# which perfectly matches a 50 bp region in the 3'UTR of XPA
# Presumably, the XPA cDNA sequence used to complement XPA in XP20S-XPA does
# not have the full 3'UTR.
probe <- "ILMN_1787591";

e <- exprs(x);
pheno <- pData(x);

pheno$cell <- NA;
pheno$cell[ pheno[["treatment:ch1"]] == "Human primary fibroblast GM969 transfected with shXPA" ] <- "GM969-shXPA";
pheno$cell[ pheno[["treatment:ch1"]] == "Human primary fibroblast GM969 transfected with a scrambled shXPA" ] <- "GM969-shControl";
pheno$cell[ pheno[["treatment:ch1"]] == "SV40 transformed XPA-deficient (GM04312C) cell line" ] <- "XP20S";
pheno$cell[ pheno[["treatment:ch1"]] == "Matched complemented XPA+ (GM15876A) cell line" ] <- "XP20S-XPA";

y <- data.frame(pheno, XPA = e[probe, ]);

ggplot(y, aes(x=cell, y=XPA)) +
	geom_jitter(stat="identity", width=0.1) + theme_bw()

