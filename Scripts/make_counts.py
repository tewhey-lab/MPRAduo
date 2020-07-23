###Make dictionary of parsed barcode - oligo pairs
###go through fastq file to pull barcodes and oligos for each


import pandas as pd
import os
import sys
import re
import time
import pathlib
import gzip
argv = sys.argv

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

fastqfile = argv[1]
dictfile = argv[2]
out_id = argv[3]
rev = int(argv[4])

current_path = os.getcwd()

BC_list = pd.read_table(dictfile, header=None, index_col=False, usecols=[0,1])
BC_dict = dict(zip(BC_list[0], BC_list[1]))

with open("%s/%s.match" % (current_path, out_id), "w") as match_oligo:
    with open("%s/%s.reject.fastq" % (current_path, out_id), "w") as reject_fastq:
        with gzip.open(fastqfile, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                seq_only = record.seq
                if rev == 1:
                    seq_only = seq_only.reverse_complement()
                    seq_only = str(seq_only)
                    bc_seq = seq_only[1:21]
                else:
                    seq_only = str(seq_only)
                    bc_seq = seq_only[0:20]

                # Check for barcode presence in the dictionary
                if bc_seq in BC_dict:
                    match_oligo.write("%s\t%s\t%s\n" % (record.name,bc_seq,BC_dict[bc_seq]))
                if bc_seq not in BC_dict:
                    SeqIO.write(record, reject_fastq, "fastq")
