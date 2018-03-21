SAMPLE=$1
OUTDIR=$2
NTHREADS=$3
WT_REF=./genomes/cxcr4/cxcr4.fa
EDITED_REF=./genomes/cxcr4/cxcr4-edited.fa
DONORSEQ=GGATCCAAGCTT

# Activate the correct python virtualenv
#source /share/pipeline/anaconda/bin/activate

mkdir -p $OUTDIR
python crispr-split-fq.py ${SAMPLE}_R2.fastq.gz $DONORSEQ $OUTDIR

bwa mem -t $NTHREADS -M $WT_REF ${OUTDIR}/${SAMPLE}_R2-wt.fastq.gz | \
samtools sort -@ NTHREADS -O BAM --reference $WT_REF > ${OUTDIR}/${SAMPLE}-wt.bam && \
samtools index -b ${OUTDIR}/${SAMPLE}-wt.bam && \
freebayes -f $WT_REF ${OUTDIR}/${SAMPLE}-wt.bam -= > ${OUTDIR}/${SAMPLE}-wt.vcf

bwa mem -t $NTHREADS -M $EDITED_REF ${OUTDIR}/${SAMPLE}_R2-edited.fastq.gz | \
samtools sort -@ $NTHREADS -O BAM --reference $EDITED_REF > ${OUTDIR}/${SAMPLE}-edited.bam && \
samtools index -b ${OUTDIR}/${SAMPLE}-edited.bam && \
freebayes -f $EDITED_REF ${OUTDIR}/${SAMPLE}-edited.bam -= > ${OUTDIR}/${SAMPLE}-edited.vcf
