BATCHNAME=$1
INDIR=$2

CURDIR=$PWD
cd $INDIR


echo "Batch\tSample\tReferenceReads\tEditedReads\tReferenceVariants\tEditedVariants";

for sample in $(ls *-wt.bam | cut -d"-" -f1); do
    ref_reads=$(samtools idxstats ${sample}-wt.bam | head -n1 | cut -f3);
    edited_reads=$(samtools idxstats ${sample}-edited.bam | head -n1 | cut -f3);

    ref_vars=$(vcffilter ${sample}-wt.vcf -f "AC > 0" -f "QUAL > 30" | grep -v "^#" | wc -l);
    edited_vars=$(vcffilter ${sample}-edited.vcf -f "AC > 0" -f "QUAL > 30" | grep -v "^#" | wc -l);
    
    echo "${BATCHNAME}\t${sample}\t${ref_reads}\t${edited_reads}\t${ref_vars}\t${edited_vars}";
done;

cd $CURDIR
