## Re-writing this to apply to single end reads instead of paired end reads

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
#input fastq file to check as first argument, input desired output directory as
#second argument
fastqfile = argv[1]
dictE = argv[2]
dictS = argv[3]
out_dir = argv[4]

start_time = time.time()

current_path = os.getcwd()

p = pathlib.Path('%s/%s' % (current_path, out_dir))
if p.exists() == False and p.is_dir() == False:
    os.mkdir(p)

#create lists of the bar codes for each library
BC_E = pd.read_table(dictE, header=None, index_col=False, usecols=[0,1])
BC_S = pd.read_table(dictS, header=None, index_col=False, usecols=[0,1])

libE_dict = dict(zip(BC_E[0], BC_E[1]))
libS_dict = dict(zip(BC_S[0], BC_S[1]))

i=0
E_found = 0
S_found = 0
count_single_no_match = 0
count_single_E_no_match = 0
count_single_S_no_match = 0
count_duo_no_match = 0
count_duo_no_match_ES = 0
count_duo_no_match_SE = 0
count_duo_ES_no_oligo = 0
count_duo_SE_no_oligo = 0
count_duo_ES_no_E = 0
count_duo_ES_no_S = 0
count_duo_SE_no_S = 0
count_duo_SE_no_E = 0
count_single = 0
count_single_E = 0
count_single_S = 0
count_duo = 0
count_duo_ES = 0
count_duo_ES_both = 0
count_duo_SE = 0
count_duo_SE_both = 0
no_oligo_match = 0
# Open write files and begin parsing the fastq file
with open("%s/%s/%s.match" % (current_path, out_dir, out_dir), "w") as fmatch:
    with open("%s/%s/%s_perfect_match.txt" % (current_path, out_dir, out_dir), "w") as perfect_oligo:
        with open("%s/%s/%s_partial_match.fastq" % (current_path, out_dir, out_dir), "w") as partial_fastq:
            with open("%s/%s/%s_partial_match.txt" % (current_path, out_dir, out_dir), "w") as partial_oligo:
                with open("%s/%s/%s_single_no_match.fastq" % (current_path, out_dir, out_dir), "w") as no_match_single:
                    with open("%s/%s/%s_duo_no_match.fastq" % (current_path, out_dir, out_dir), "w") as no_match_duo:
                        with open("%s/%s/%s_oligo_no_match.fastq" % (current_path, out_dir, out_dir), "w") as no_match_oligo:
                            with gzip.open(fastqfile, "rt") as handle:
                                for record in SeqIO.parse(handle, "fastq"):
    # Identify the read and set GFP_end to the appropriate 8 base identifier
    # If it is read 1, reverse the read
                                    seq_only = record.seq
                                    seq_only = str(seq_only)
    # print(seq_only[0:11])
                                    UMI_end = 'CATGTC'
    # Here E is associated with plasmid A, and S is associated with plasmid P
                                    ES_indic = 'CGATCTCA'
                                    SE_indic = 'GTTTCTAG'
                                    link_E = 'CGCCGA'
                                    link_S = 'CTCGGC'
    # Find the end of the UMI read in the sequence line
                                    if UMI_end in seq_only:
                                        match = seq_only.index(UMI_end)
                                    else:
                                        match = 0
                                    fill_start = match
                                    lib_E_index = 0
                                    lib_S_index = 0
    # Split the sequences based on the presence of the duo linker
    # 3 steps check for ES, then SE, then assume single if not present in either.
                                    if ES_indic in seq_only[fill_start+71:fill_start+85]:
                                        count_duo += 1
                                        count_duo_ES += 1
                                        if link_E in seq_only[45:67]:
                                            lib_E_index = seq_only.index(link_E)
                                            E_found += 1
                                        else:
                                            lib_E_index = -1
                                        if link_S in seq_only[93:105]:
                                            lib_S_index = seq_only.index(link_S)
                                            S_found += 1
                                        else:
                                            lib_S_index = -1
                                        if lib_E_index in range(45,67) and lib_S_index in range(93,105):
                                            k = seq_only[lib_E_index+6:lib_E_index+16]
                                            g = seq_only[lib_S_index+6:lib_S_index+26]
                                            bc_E = Seq(k)
                                            bc_S = Seq(g)
                                            bc_E = bc_E.reverse_complement()
                                            bc_S = bc_S.reverse_complement()
                                            bc_E = str(bc_E)
                                            bc_S = str(bc_S)
                                            if bc_E in libE_dict and bc_S in libS_dict:
                                                count_duo_ES_both += 1
                                                perfect_oligo.write("%s^%s\t%s^%s\n" % (bc_E,bc_S,libE_dict[bc_E],libS_dict[bc_S]))
                                                fmatch.write("%s\t%s^%s\t%s^%s\tES\n" % (record.name,bc_E,bc_S,libE_dict[bc_E],libS_dict[bc_S]))
                                            elif bc_E in libE_dict and bc_S not in libS_dict:
                                                count_duo_ES_no_S += 1
                                                partial_oligo.write("%s\t%s\t%s\t-\n" % (bc_E,bc_S,libE_dict[bc_E]))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                            elif bc_E not in libE_dict and bc_S in libS_dict:
                                                count_duo_ES_no_E += 1
                                                partial_oligo.write("%s\t%s\t-\t%s\n" % (bc_E,bc_S,libS_dict[bc_S]))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                            elif bc_E not in libE_dict and bc_S not in libS_dict:
                                                no_oligo_match += 1
                                                count_duo_ES_no_oligo += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        elif lib_E_index in range(45,62) and lib_S_index not in range(93,105):
                                            count_duo_no_match_ES += 1
                                            count_duo_ES_no_S += 1
                                        elif lib_E_index not in range(45,62) and lib_S_index in range(93,105):
                                            count_duo_no_match_ES += 1
                                            count_duo_ES_no_E += 1
                                        else:
                                            count_duo_no_match_ES += 1
                                    elif SE_indic in seq_only[fill_start+81:fill_start+95]:
                                        count_duo += 1
                                        count_duo_SE += 1
                                        if link_S in seq_only[45:62]:
                                            lib_S_index = seq_only.index(link_S)
                                        else:
                                            lib_S_index = -1
                                        if link_E in seq_only[93:105]:
                                            lib_E_index = seq_only.index(link_E)
                                        else:
                                            lib_E_index = -1
                                        if lib_S_index in range(45,62) and lib_E_index in range(93,105):
                                            k = seq_only[lib_E_index+6:lib_E_index+16]
                                            g = seq_only[lib_S_index+6:lib_S_index+26]
                                            bc_E = Seq(k)
                                            bc_S = Seq(g)
                                            bc_E = bc_E.reverse_complement()
                                            bc_S = bc_S.reverse_complement()
                                            bc_E = str(bc_E)
                                            bc_S = str(bc_S)
                                            if bc_E in libE_dict and bc_S in libS_dict:
                                                count_duo_ES_both += 1
                                                perfect_oligo.write("%s^%s\t%s^%s\n" % (bc_S,bc_E,libS_dict[bc_S],libE_dict[bc_E]))
                                                fmatch.write("%s\t%s^%s\t%s^%s\tSE\n" % (record.name,bc_S,bc_E,libS_dict[bc_S],libE_dict[bc_E]))
                                            elif bc_S in libS_dict and bc_E not in libE_dict:
                                                count_duo_SE_no_E += 1
                                                partial_oligo.write("%s\t%s\t%s\t-\n" % (bc_S,bc_E,libS_dict[bc_S]))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                            elif bc_S not in libS_dict and bc_E in libE_dict:
                                                count_duo_SE_no_S += 1
                                                partial_oligo.write("%s\t%s\t-\t%s\n" % (bc_S,bc_E,libE_dict[bc_E]))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                            elif bc_E not in libE_dict and bc_S not in libS_dict:
                                                no_oligo_match += 1
                                                count_duo_SE_no_oligo += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        elif lib_S_index in range(45,62) and lib_E_index not in range(93,105):
                                            count_duo_no_match_SE += 1
                                            count_duo_SE_no_S += 1
                                        elif lib_S_index not in range(45,62) and lib_E_index in range(93,105):
                                            count_duo_no_match_ES += 1
                                            count_duo_SE_no_E += 1
                                        else:
                                            count_duo_no_match_SE += 1
                                    else:
                                        count_single += 1
                                        if link_E in seq_only[fill_start+38:fill_start+50]:
                                            count_single_E += 1
                                            lib_E_index = seq_only.index(link_E)
                                            k = seq_only[lib_E_index+6:lib_E_index+16]
                                            bc_E = Seq(k)
                                            bc_E = bc_E.reverse_complement()
                                            bc_E = str(bc_E)
                                            if bc_E in libE_dict:
                                                perfect_oligo.write("%s\t%s\n" % (bc_E, libE_dict[bc_E]))
                                                fmatch.write("%s\t%s\t%s\tE\n" % (record.name,bc_E,libE_dict[bc_E]))
                                            else:
                                                no_oligo_match += 1
                                                count_single_E_no_match += 1
                                                count_single_no_match += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        elif link_S in seq_only[fill_start+38:fill_start+50]:
                                            lib_S_index = seq_only.index(link_S)
                                            g = seq_only[lib_S_index+6:lib_S_index+26]
                                            bc_S = Seq(g)
                                            bc_S = bc_S.reverse_complement()
                                            bc_S = str(bc_S)
                                            if bc_S in libS_dict:
                                                perfect_oligo.write("%s\t%s\n" % (bc_S, libS_dict[bc_S]))
                                                fmatch.write("%s\t%s\t%s\tS\n" % (record.name,bc_S,libS_dict[bc_S]))
                                            else:
                                                no_oligo_match += 1
                                                count_single_S_no_match += 1
                                                count_single_no_match += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                    if i%100 == 0:
                                        print(i)
                                    i += 1

print(i)
print(E_found)
print(S_found)
if count_duo_SE == 0:
    count_duo_SE = -1
if count_duo_ES == 0:
    count_duo_ES = -1
if count_single_S == 0:
    count_single_S = -1
if count_single_E == 0:
    count_single_E = -1
print("%s seconds" % (time.time() - start_time))
with open ("%s/%s/statistics.txt"  % (current_path, out_dir), "w") as run_stats:
    run_stats.write("Total Records Checked: %f \n" % float(i))
    run_stats.write("Percent Matched to Library: %f \n" % float(1-((count_duo_no_match+count_duo_ES_no_S+count_duo_ES_no_E+count_duo_SE_no_E+count_duo_SE_no_S+count_single_no_match)/i)))
    run_stats.write("Percent Duo Matched: %f \n" % float(1-((count_duo_no_match+count_duo_ES_no_S+count_duo_ES_no_E+count_duo_SE_no_E+count_duo_SE_no_S)/count_duo)))
    run_stats.write("Percent Single Matched: %f \n" % float(1-(count_single_no_match/count_single)))
    run_stats.write("Percent Matched with No Oligos: %f\n" % float(no_oligo_match/((count_duo-count_duo_no_match)+(count_single-count_single_no_match))))
    run_stats.write("Percent ES Mis-Match: %f\t%f\t%f\n" % (float(count_duo_no_match_ES/count_duo_ES), float(count_duo_no_match_ES), float(count_duo_ES)))
    run_stats.write("Percent SE Mis-Match: %f\t%f\t%f\n" % (float(count_duo_no_match_SE/count_duo_SE), float(count_duo_no_match_SE), float(count_duo_SE)))
    run_stats.write("Single E: \n \tPercent Oligos Matched: %f\n \tPercent Oligos Not Matched: %f\n" % (float(1-(count_single_E_no_match/count_single_E)), float(count_single_E_no_match/count_single_E)))
    run_stats.write("Single S: \n \tPercent Oligos Matched: %f\n \tPercent Oligos Not Matched: %f\n" % (float(1-(count_single_S_no_match/count_single_S)), float(count_single_S_no_match/count_single_S)))
    run_stats.write("ES: \n \tPercent Both Oligos Matched: %f\n \tPercent No Oligos Matched: %f\n \tPercent only E oligo matched: %f\n \tPercent only S oligo matched: %f\n" % (float(count_duo_ES_both/count_duo_ES), float(count_duo_ES_no_oligo/count_duo_ES), float(count_duo_ES_no_S/count_duo_ES), float(count_duo_ES_no_E/count_duo_ES)))
    run_stats.write("SE: \n \tPercent Both Oligos Matched: %f\n \tPercent No Oligos Matched: %f\n \tPercent only E oligo matched: %f\n \tPercent only S oligo matched: %f\n" % (float(count_duo_SE_both/count_duo_SE), float(count_duo_SE_no_oligo/count_duo_SE), float(count_duo_SE_no_S/count_duo_SE), float(count_duo_SE_no_E/count_duo_SE)))
