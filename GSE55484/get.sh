accession=GSE55484

url=ftp://ftp.ncbi.nlm.nih.gov/geo/series
outdir=raw

mkdir -p $outdir
cd $outdir

# download and extract main raw archive file
group=${accession%???}nnn
file=${accession}_RAW.tar
curl -o ${file} ${url}/${group}/${accession}/suppl/$file
tar -xf ${file} && rm ${file}
gunzip *.gz

# extra file
extra_file=${accession}_non-normalized
extra_gz=${extra_file}.txt.gz
curl -o ${extra_gz} ${url}/${group}/${accession}/suppl/${extra_gz}
gunzip ${extra_gz}

# custom preprocessing
# remove extra trailing whitespace
# remove comments (third line is malformed
sed 's/\t\t\+//' ${extra_file}.txt |
	sed '1,3d' > ${extra_file}.tsv

cd -
