wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55484/suppl/GSE55484_non-normalized.txt.gz

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55484/suppl/GSE55484_RAW.tar
tar -xf GSE55484_RAW.tar
rm GSE55484_RAW.tar

gunzip *.txt.gz
sed -i 's/\t\t\+//' GSE55484_non-normalized.txt
