#!/usr/bin/env python
"""Out put how many sequences are found in common between 
two sam files that only contain proper match sequences.
The iput to this script should be two files on the command line.
These two files should just be the alignment section of a SAM file that contains
only the propper pairs. To get one of these files from a BAM files:
samtools view -f 2 M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.mp.lutea_align.bam > proper_pairs.bwa.txt
"""

import sys

class CompareSAMS:
    def __init__(self):
        # get list of seq names in file one
        with open(sys.argv[-2], 'r') as f:
            self.seq_names_set_one = set([line.split()[0] for line in f])
        # get list of seq names in file two
        with open(sys.argv[-1], 'r') as f:
            self.seq_names_set_two = set([line.split()[0] for line in f])
        self.len_one = len(self.seq_names_set_one)
        self.len_two = len(self.seq_names_set_two)
        self.seqs_in_common = len(self.seq_names_set_one.intersection(self.seq_names_set_two))
        self.unique_seqs_one = self.len_one - self.seqs_in_common
        self.unique_seqs_two = self.len_two - self.seqs_in_common
        self.total_unique_seqs = len(self.seq_names_set_one.union(self.seq_names_set_two))
        # print stats
        print(
            f'{self.len_one} proper pairs in {sys.argv[-2]}\n'
            f'{self.len_two} proper pairs in {sys.argv[-1]}\n'
            f'{self.len_one-self.seqs_in_common} unique proper pairs in {sys.argv[-2]}\n'
            f'{self.len_two-self.seqs_in_common} unique proper pairs in {sys.argv[-1]}\n'
            f'{self.seqs_in_common} proper pairs found in common ({self.seqs_in_common/self.total_unique_seqs})\n'
            )
cs = CompareSAMS()