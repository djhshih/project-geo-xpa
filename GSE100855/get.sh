accession=GSE100855

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

cd -
