from collections import namedtuple
import csv
import gzip
from itertools import izip
from path import Path
import sys


FQ = namedtuple('FQ', ['header', 'seq', 'dummy', 'qual'])


def fqiter(fname):
    with gzip.open(fname) as fp:
        # Read entire file in memory. Faster than readline line at a time over NFS!
        data = iter(fp.readlines())
        for header in data:
            seq = data.next()
            dummy = data.next()
            qual = data.next()
            yield FQ(header=header, seq=seq, dummy=dummy, qual=qual)


def fq_add_modifier(fname, mod):
    if fname.endswith('fastq.gz'):
        return "%s-%s.fastq.gz" % (fname[:-9], mod)
    elif fname.endswith('fastq'):
        return "%s-%s.fastq" % (fname[:-6], mod)
    raise ValueError("Unexpected file name: %s" % fname)


def revcomp(seq):
    rev = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return ''.join(rev[s] for s in reversed(seq))


def split_fq(fname, donor_seq, outdir):
    outdir = Path(outdir)
    donor_seq_rc = revcomp(donor_seq)
    wt_fname = outdir / fq_add_modifier(fname, 'wt')
    edited_fname = outdir / fq_add_modifier(fname, 'edited')

    with gzip.open(wt_fname, 'wb') as wt_fq, gzip.open(edited_fname, 'wb') as edited_fq:
        for read in fqiter(fname):
            edited_seq = donor_seq in read.seq or donor_seq_rc in read.seq
            out_fq = edited_fq if edited_seq else wt_fq
            out_fq.writelines((read.header, read.seq, read.dummy, read.qual))

    return wt_fname, edited_fname


if __name__ == '__main__':
    split_fq(*sys.argv[1:])
