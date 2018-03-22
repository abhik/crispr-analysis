"""
Split a FASTQ file into multiple files based on inline barcodes on both reads.

Usage:
    inline-demux.py <READ1_FASTQ> <READ2_FASTQ> <BARCODE_LIST> <OUTPUT_DIRECTORY> <OUTPUT_PREFIX>

READ1_FASTQ and READ2_FASTQ: 
    Paths to the input fastq files

BARCODE_LIST: 
    A tab-delimited file with following columns: sample name, read1 barcode, read2 barcode
    file should not have a header row.

OUTPUT_DIRECTORY:
    Directory where output files should be created.

OUTPUT_PREFIX:
    A prefix to add to all output files.
    Output files will be named like <OUTPUT_PREFIX>_<SAMPLE_NAME>_R1.fastq.gz
    
---------------------------------------------------------------------------------------------------------
Given fastq files named TUMOR1_S2_R1.fastq.gz and TUMOR1_S2_R2.fastq.gz, this is an example usage:

> python /share/bli/inline-demux.py TUMOR1_S2_R1.fastq.gz TUMOR1_S2_R2.fastq.gz barcodes.txt demuxed_data TUMOR1
"""


from collections import namedtuple
import csv
import gzip
from itertools import izip
from path import Path
import sys


FQ = namedtuple('FQ', ['header', 'seq', 'dummy', 'qual'])


def fqiter(fname):
    with gzip.open(fname) as fp:
	# Read entire file in memory. 
    # Faster than reading line at a time over NFS. Won't be true when reading files locally.
	data = iter(fp.readlines())
    for header in data:
        seq = data.next()
        dummy = data.next()
        qual = data.next()
        yield FQ(header=header, seq=seq, dummy=dummy, qual=qual)


def setup_outfqs(bcinfo, outdir, outfq_prefix):
    outdir = Path(outdir)
    sample2outfqs =  {
        sample: (gzip.open(outdir / "%s_%s_R1.fastq.gz" % (outfq_prefix, sample), 'wb'),
                 gzip.open(outdir / "%s_%s_R2.fastq.gz" % (outfq_prefix, sample), 'wb'))
        for sample, bc1, bc2 in bcinfo}

    sample2outfqs[None] = (gzip.open(outdir / "%s_Undetermined_R1.fastq.gz" % outfq_prefix, 'wb'),
                           gzip.open(outdir / "%s_Undetermined_R2.fastq.gz" % outfq_prefix, 'wb'))
    return sample2outfqs


def demux(fname1, fname2, bcinfo_fname, outdir, outfq_prefix):
    with open(bcinfo_fname) as fp:
        bcinfo = list(csv.reader(fp, delimiter="\t"))

    bc2sample = {(bc1, bc2): sample for sample, bc1, bc2 in bcinfo}
    sample2outfqs = setup_outfqs(bcinfo, outdir, outfq_prefix)
    sample_counts = {sample: 0 for sample in sample2outfqs.keys()}
    bclen = len(bcinfo[0][1])

    # demux reads with perfect barcode matches
    for fq1, fq2 in izip(fqiter(fname1), fqiter(fname2)):
        bc1 = fq1.seq[:bclen]
        bc2 = fq2.seq[:bclen]
        sample = bc2sample.get((bc1, bc2))
        sample_counts[sample] += 1
        outfq1, outfq2 = sample2outfqs.get(sample)
        if outfq1:
            for infq, outfq in zip((fq1, fq2), (outfq1, outfq2)):
                outfq.write(infq.header)
                outfq.write(infq.seq[bclen:])
                outfq.write(infq.dummy)
                outfq.write(infq.qual[bclen:])

    # close outfqs
    for outfq1, outf12 in sample2outfqs.itervalues():
        outfq1.close()
        outfq2.close()

    # print counts
    with open(Path(outdir) / "%s.demux-stats.txt" % outfq_prefix, 'w') as fp:
        for sample, count in sample_counts.iteritems():
            fp.write("%s\t%d\n" % (sample, count))

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print __doc__
        sys.exit()
    demux(*sys.argv[1:])
