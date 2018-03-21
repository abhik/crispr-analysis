Scripts to compute editing efficiency of single-cell CRISPR experiments.

## Example workflow for analyzing CRISPR data

Given fastq files in `./fastqs` and barcodes in `./6bp-barcodes.txt`, this set of commands will analyze HDR and NHEJ edits.

```sh
mkdir crispr && cp ./fastqs/AM_*.fastq.gz ./crispr
cp 6bp-barcodes.txt ./crispr

# Demux using inline barcodes (can skip this step if fastqs already demuxed)
cd crispr
ls *R1*.gz | cut -d"_" -f1,2 > samples.txt
mkdir demuxed
cat samples.txt | parallel python ../inline-demux.py {}_*_R1_001.fastq.gz {}_*_R2_001.fastq.gz 6bp-barcodes.txt demuxed {}

# Generate crispr results
mkdir -p alignments && cat samples.txt | parallel sh -x ../crispr.sh {} alignments 2
cd alignments && sh ../crispr-collect-results.sh 161216 . > crispr-results.txt
```
