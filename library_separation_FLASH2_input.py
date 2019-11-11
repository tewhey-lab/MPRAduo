import pandas as pd
import os
import sys
import re
import time
import pathlib
argv = sys.argv
sys.path.append("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/")

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
#input fastq file to check as first argument, input desired output directory as
#second argument
fastqfile = argv[1]
out_dir = argv[2]

start_time = time.time()

p = pathlib.Path('/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s' % out_dir)
if p.exists() == False and p.is_dir() == False:
    os.mkdir(p)

A_tag_asso = '/projects/tewhey-lab/mourik/MPRAduo_2nd_assembly/MPRAduo_LibA_pilot_20190415.merged.match.enh.mapped.barcode.ct'
P_tag_asso = '/projects/tewhey-lab/mourik/MPRAduo_2nd_assembly/MPRAduo_LibP_pilot.merged.match.enh.mapped.barcode.ct.parsed'

#create lists of the bar codes for each library
BC_A = pd.read_table(A_tag_asso, header=None, index_col=False, usecols=[0,1])
BC_P = pd.read_table(P_tag_asso, header=None, index_col=False, usecols=[0,1])

libA_dict = dict(zip(BC_A[0], BC_A[1]))
libP_dict = dict(zip(BC_P[0], BC_P[1]))

i=0
count_single_no_match = 0
count_single_A_no_match = 0
count_single_P_no_match = 0
count_duo_no_match = 0
count_duo_no_match_AiP = 0
count_duo_no_match_PiA = 0
count_duo_AiP_no_oligo = 0
count_duo_PiA_no_oligo = 0
count_duo_AiP_no_A = 0
count_duo_AiP_no_P = 0
count_duo_PiA_no_P = 0
count_duo_PiA_no_A = 0
count_single = 0
count_single_A = 0
count_single_P = 0
count_duo = 0
count_duo_AiP = 0
count_duo_AiP_both = 0
count_duo_PiA = 0
count_duo_PiA_both = 0
no_oligo_match = 0
# Open write files and begin parsing the fastq file
with open("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/perfect_match.txt" % out_dir, "w") as perfect_oligo:
    with open("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/partial_match.fastq" % out_dir, "w") as partial_fastq:
        with open("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/partial_match.txt" % out_dir, "w") as partial_oligo:
            with open("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/single_no_match.fastq" % out_dir, "w") as no_match_single:
                with open("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/duo_no_match.fastq" % out_dir, "w") as no_match_duo:
                    with open("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/oligo_no_match.fastq" % out_dir, "w") as no_match_oligo:
                        with open(fastqfile, "rt") as handle:
                            for record in SeqIO.parse(handle, "fastq"):
# Identify the read and set GFP_end to the appropriate 8 base identifier
# If it is read 1, reverse the read
                                seq_only = record.seq
                                seq_only = seq_only.reverse_complement()
                                seq_only = str(seq_only)
# print(seq_only[0:11])
                                GFP_end = 'AATAATA'
                                link_A = 'TCTAGA'
                                link_P = 'CTGACT'
# Find the end of the GFP read in the sequence line
                                if GFP_end in seq_only:
                                    match = seq_only.index(GFP_end)
                                else:
                                    match = 0
                                fill_start = match
# Split the sequences up by the lenght in relation to the GFP_end, less than 110
# bases to the end of the sequence is a single
                                if len(seq_only[fill_start:]) < 110:
                                    lib_A_index = 0
                                    lib_P_index = 0
                                    count_single += 1
                                    if link_A in seq_only[fill_start+19:fill_start+28]:
                                        lib_A_index = seq_only.index(link_A)
                                    if link_P in seq_only[fill_start+19:fill_start+28]:
                                        lib_P_index = seq_only.index(link_P)
                                    if lib_A_index > 0 and lib_P_index == 0:
                                        count_single_A += 1
                                        g = seq_only[lib_A_index+6:lib_A_index+26]
                                        if g in libA_dict:
                                            perfect_oligo.write("%s\t%s\n" % (g,libA_dict[g]))
                                        if g not in libA_dict:
                                            # matched_all_oligo.write("%s\n" % g)
                                            no_oligo_match += 1
                                            count_single_A_no_match += 1
                                            SeqIO.write(record, no_match_oligo, "fastq")
                                    if lib_P_index > 0 and lib_A_index == 0:
                                        count_single_P += 1
                                        k = seq_only[lib_P_index+6:lib_P_index+16]
                                        if k in libP_dict:
                                            perfect_oligo.write("%s\t%s\n" % (k, libP_dict[k]))
                                        if k not in libP_dict:
                                            # matched_all_oligo.write("%s\n" % k)
                                            no_oligo_match += 1
                                            count_single_P_no_match += 1
                                            SeqIO.write(record, no_match_oligo, "fastq")
                                    if lib_A_index == 0 and lib_P_index == 0:
                                        count_single_no_match += 1
                                        SeqIO.write(record, no_match_single, "fastq")
                                if len(seq_only[fill_start:]) >= 110:
                                    lib_A_index = 0
                                    lib_P_index = 0
                                    lib_AiP_index = 0
                                    lib_PiA_index = 0
                                    count_duo += 1
                                    if link_A in seq_only[fill_start+19:fill_start+28]:
                                        lib_A_index = seq_only.index(link_A)
                                    if link_P in seq_only[lib_A_index+55:lib_A_index+64]:
                                        seq_only_AiP = seq_only[lib_A_index+55:]
                                        lib_AiP_index = seq_only_AiP.index(link_P)
                                    if link_P in seq_only[fill_start+19:fill_start+28]:
                                        lib_P_index = seq_only.index(link_P)
                                    if link_A in seq_only[lib_P_index+55:lib_P_index+64]:
                                        seq_only_PiA = seq_only[lib_P_index+55:]
                                        lib_PiA_index = seq_only_PiA.index(link_A)
                                    if lib_A_index > 0 and lib_AiP_index > 0:
                                        count_duo_AiP += 1
                                        g = seq_only[lib_A_index+6:lib_A_index+26]
                                        k = seq_only_AiP[lib_AiP_index+6:lib_AiP_index+16]
                                        if g in libA_dict and k in libP_dict:
                                            count_duo_AiP_both += 1
                                            perfect_oligo.write("%s^%s\t%s^%s\n" % (g,k,libA_dict[g],libP_dict[k]))
                                        if g in libA_dict and k not in libP_dict:
                                            partial_oligo.write("%s\t%s\t%s\t-\n" % (g,k,libA_dict[g]))
                                            SeqIO.write(record, partial_fastq, "fastq")
                                            count_duo_AiP_no_P += 1
                                            # matched_all_oligo.write("%s^%s\t%s^\n" % (g,k,libA_dict[g]))
                                        if g not in libA_dict and k in libP_dict:
                                            partial_oligo.write("%s\t%s\t-\t%s\n" % (g,k,libP_dict[k]))
                                            SeqIO.write(record, partial_fastq, "fastq")
                                            # matched_all_oligo.write("%s^%s\t^%s\n" % (g,k,libP_dict[k]))
                                            count_duo_AiP_no_A += 1
                                        if g not in libA_dict and k not in libP_dict:
                                            # matched_all_oligo.write("%s^%s\t^\n" % (g,k))
                                            no_oligo_match += 1
                                            count_duo_AiP_no_oligo += 1
                                            SeqIO.write(record, no_match_oligo, "fastq")
                                    if lib_P_index > 0 and lib_PiA_index > 0:
                                        count_duo_PiA += 1
                                        g = seq_only[lib_P_index+6:lib_P_index+16]
                                        k = seq_only_PiA[lib_PiA_index+6:lib_PiA_index+26]
                                        if g in libP_dict and k in libA_dict:
                                            count_duo_PiA_both += 1
                                            perfect_oligo.write("%s^%s\t%s^%s\n" % (g,k,libP_dict[g],libA_dict[k]))
                                        if g in libP_dict and k not in libA_dict:
                                            partial_oligo.write("%s\t%s\t%s\t-\n" % (g,k,libA_dict[g]))
                                            SeqIO.write(record, partial_fastq, "fastq")
                                            count_duo_PiA_no_A += 1
                                            # matched_all_oligo.write("%s^%s\t%s^\n" % (g,k,libP_dict[g]))
                                        if g not in libP_dict and k in libA_dict:
                                            partial_oligo.write("%s\t%s\t-\t%s\n" % (g,k,libP_dict[k]))
                                            SeqIO.write(record, partial_fastq, "fastq")
                                            count_duo_PiA_no_P += 1
                                            # matched_all_oligo.write("%s^%s\t^%s\n" % (g,k,libA_dict[k]))
                                        if g not in libP_dict and k not in libA_dict:
                                            # matched_all_oligo.write("%s^%s\t^\n" % (g,k))
                                            no_oligo_match += 1
                                            count_duo_PiA_no_oligo += 1
                                            SeqIO.write(record, no_match_oligo, "fastq")
                                    if lib_A_index == 0 and lib_AiP_index > 0 or lib_A_index > 0 and lib_AiP_index == 0:
                                        count_duo_no_match += 1
                                        count_duo_no_match_AiP += 1
                                        SeqIO.write(record, no_match_duo, "fastq")
                                    if lib_P_index == 0 and lib_PiA_index > 0 or lib_P_index > 0 and lib_PiA_index == 0:
                                        count_duo_no_match += 1
                                        count_duo_no_match_PiA += 1
                                        SeqIO.write(record, no_match_duo, "fastq")
                                    if lib_P_index == 0 and lib_A_index == 0 and lib_PiA_index == 0 and lib_AiP_index == 0:
                                        count_duo_no_match += 1
                                        SeqIO.write(record, no_match_duo, "fastq")
                                if i%100 == 0:
                                    print(i)
                                i += 1
                                # if i > 550:
                                #     break

print("%s seconds" % (time.time() - start_time))
with open ("/projects/tewhey-lab/mourik/MPRAduo_Hannah/output/%s/statistics.txt"  % out_dir, "w") as run_stats:
    run_stats.write("Total Records Checked: %f \n" % float(i))
    run_stats.write("Percent Matched to Library: %f \n" % float(1-((count_duo_no_match+count_single_no_match)/i)))
    run_stats.write("Percent Duo Matched: %f \n" % float(1-(count_duo_no_match/count_duo)))
    run_stats.write("Percent Single Matched: %f \n" % float(1-(count_single_no_match/count_single)))
    run_stats.write("Percent Matched with No Oligos: %f\n" % float(no_oligo_match/((count_duo-count_duo_no_match)+(count_single-count_single_no_match))))
    run_stats.write("Percent A into P Mis-Match: %f\n" % float(count_duo_no_match_AiP/count_duo_AiP))
    run_stats.write("Percent P into A Mis-Match: %f\n" % float(count_duo_no_match_PiA/count_duo_PiA))
    run_stats.write("Single A: \n \tPercent Oligos Matched: %f\n \tPercent Oligos Not Matched: %f\n" % (float(1-(count_single_A_no_match/count_single_A)), float(count_single_A_no_match/count_single_A)))
    run_stats.write("Single P: \n \tPercent Oligos Matched: %f\n \tPercent Oligos Not Matched: %f\n" % (float(1-(count_single_P_no_match/count_single_P)), float(count_single_P_no_match/count_single_P)))
    run_stats.write("A into P: \n \tPercent Both Oligos Matched: %f\n \tPercent No Oligos Matched: %f\n \tPercent only A oligo matched: %f\n \tPercent only P oligo matched: %f\n" % (float(count_duo_AiP_both/count_duo_AiP), float(count_duo_AiP_no_oligo/count_duo_AiP), float(count_duo_AiP_no_P/count_duo_AiP), float(count_duo_AiP_no_A/count_duo_AiP)))
    run_stats.write("P into A: \n \tPercent Both Oligos Matched: %f\n \tPercent No Oligos Matched: %f\n \tPercent only A oligo matched: %f\n \tPercent only P oligo matched: %f\n" % (float(count_duo_PiA_both/count_duo_PiA), float(count_duo_PiA_no_oligo/count_duo_PiA), float(count_duo_PiA_no_P/count_duo_PiA), float(count_duo_PiA_no_A/count_duo_PiA)))
