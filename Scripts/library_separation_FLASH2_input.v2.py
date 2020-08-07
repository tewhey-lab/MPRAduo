## This should be used with paired end reads

import pandas as pd
import os
import sys
import re
import time
import pathlib
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
count_single_no_match = 0
count_single_E_no_match = 0
count_single_S_no_match = 0
count_duo_no_match = 0
count_duo_no_match_SE = 0
count_duo_no_match_ES = 0
count_duo_SE_no_oligo = 0
count_duo_ES_no_oligo = 0
count_duo_SE_no_E = 0
count_duo_SE_no_S = 0
count_duo_ES_no_S = 0
count_duo_ES_no_E = 0
count_single = 0
count_single_E = 0
count_single_S = 0
count_duo = 0
count_duo_SE = 0
count_duo_SE_both = 0
count_duo_ES = 0
count_duo_ES_both = 0
no_oligo_match = 0
# Open write files and begin parsing the fastq file
with open("%s.match" % out_dir, "w") as fmatch:
    with open("%s_perfect_match.txt" % out_dir, "w") as perfect_oligo:
        with open("%s_partial_match.fastq" % out_dir, "w") as partial_fastq:
            with open("%s_partial_match.txt" % out_dir, "w") as partial_oligo:
                with open("%s_single_no_match.fastq" % out_dir, "w") as no_match_single:
                    with open("%s_duo_no_match.fastq" % out_dir, "w") as no_match_duo:
                        with open("%s_oligo_no_match.fastq" % out_dir, "w") as no_match_oligo:
                            with open(fastqfile, "rt") as handle:
                                for record in SeqIO.parse(handle, "fastq"):
    # Identify the read and set GFP_end to the appropriate 8 base identifier
    # If it is read 1, reverse the read
                                    seq_only = record.seq
                                    seq_only = seq_only.reverse_complement()
                                    seq_only = str(seq_only)
    # print(seq_only[0:11])
                                    GFP_end = 'AATAATA'
                                    link_S = 'TCTAGA'
                                    link_E = 'CTGACT'
    # Find the end of the GFP read in the sequence line
                                    if GFP_end in seq_only:
                                        match = seq_only.index(GFP_end)
                                    else:
                                        match = 0
                                    fill_start = match
    # Split the sequences up by the lenght in relation to the GFP_end, less than 110
    # bases to the end of the sequence is a single
                                    if len(seq_only[fill_start:]) < 110:
                                        lib_E_index = 0
                                        lib_S_index = 0
                                        count_single += 1
                                        if link_E in seq_only[fill_start+19:fill_start+28]:
                                            lib_E_index = seq_only.index(link_E)
                                        if link_S in seq_only[fill_start+19:fill_start+28]:
                                            lib_S_index = seq_only.index(link_S)
                                        if lib_S_index > 0 and lib_E_index == 0:
                                            count_single_S += 1
                                            g = seq_only[lib_S_index+6:lib_S_index+26]
                                            if g in libS_dict:
                                                perfect_oligo.write("%s\t%s\n" % (g,libS_dict[g]))
                                                fmatch.write("%s\t%s\t%s\tS\n" % (record.name,g,libS_dict[g]))
                                            if g not in libS_dict:
                                                # matched_all_oligo.write("%s\n" % g)
                                                no_oligo_match += 1
                                                count_single_S_no_match += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        if lib_E_index > 0 and lib_S_index == 0:
                                            count_single_E += 1
                                            k = seq_only[lib_E_index+6:lib_E_index+16]
                                            if k in libE_dict:
                                                perfect_oligo.write("%s\t%s\n" % (k, libE_dict[k]))
                                                fmatch.write("%s\t%s\t%s\tE\n" % (record.name,k,libE_dict[k]))
                                            if k not in libE_dict:
                                                # matched_all_oligo.write("%s\n" % k)
                                                no_oligo_match += 1
                                                count_single_E_no_match += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        if lib_S_index == 0 and lib_E_index == 0:
                                            count_single_no_match += 1
                                            SeqIO.write(record, no_match_single, "fastq")
                                    if len(seq_only[fill_start:]) >= 110:
                                        lib_E_index = 0
                                        check_E_index = 0
                                        check_SE_index = 0
                                        lib_S_index = 0
                                        check_S_index = 0
                                        check_ES_index = 0
                                        lib_SE_index = 0
                                        lib_ES_index = 0
                                        count_duo += 1
                                        if link_S in seq_only[fill_start+19:fill_start+31]:
                                            lib_S_index = seq_only.index(link_S)
                                            check_S_index = 1
                                            # count_duo_AiP += 1
                                        if link_E in seq_only[fill_start+19:fill_start+31]:
                                            lib_E_index = seq_only.index(link_E)
                                            check_E_index = 1
                                        if link_E in seq_only[lib_S_index+55:lib_S_index+67]:
                                            seq_only_ES = seq_only[lib_S_index+55:]
                                            lib_ES_index = seq_only_ES.index(link_E)
                                            check_ES_index = 1
                                        if link_S in seq_only[lib_E_index+45:lib_E_index+57]:
                                            seq_only_SE = seq_only[lib_E_index+45:]
                                            lib_SE_index = seq_only_SE.index(link_S)
                                            check_SE_index = 1
                                        if (check_E_index == 0 and check_SE_index == 1) or (check_E_index == 1 and check_SE_index == 0):
                                            check_E_index = -1
                                        if (check_S_index == 0 and check_ES_index == 1) or (check_S_index == 1 and check_ES_index == 0):
                                            check_S_index = -1
                                        if lib_S_index > 0 and lib_ES_index > 0:
                                            count_duo_ES += 1
                                            g = seq_only[lib_S_index+6:lib_S_index+26]
                                            k = seq_only_ES[lib_ES_index+6:lib_ES_index+16]
                                            if g in libS_dict and k in libE_dict:
                                                count_duo_ES_both += 1
                                                perfect_oligo.write("%s^%s\t%s^%s\n" % (k,g,libE_dict[k],libS_dict[g]))
                                                fmatch.write("%s\t%s^%s\t%s^%s\tES\n" % (record.name,k,g,libE_dict[k],libS_dict[g]))
                                            if g in libS_dict and k not in libE_dict:
                                                partial_oligo.write("%s\t%s\t%s\t-\t%f\t%f\n" % (k,g,libS_dict[g], float(lib_S_index), float(lib_ES_index)))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                                count_duo_ES_no_E += 1
                                                # matched_all_oligo.write("%s^%s\t%s^\n" % (g,k,libA_dict[g]))
                                            if g not in libS_dict and k in libE_dict:
                                                partial_oligo.write("%s\t%s\t-\t%s\t%f\t%f\n" % (k,g,libE_dict[k], float(lib_S_index), float(lib_ES_index)))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                                # matched_all_oligo.write("%s^%s\t^%s\n" % (g,k,libP_dict[k]))
                                                count_duo_ES_no_S += 1
                                            if g not in libS_dict and k not in libE_dict:
                                                # matched_all_oligo.write("%s^%s\t^\n" % (g,k))
                                                no_oligo_match += 1
                                                count_duo_ES_no_oligo += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        if lib_E_index > 0 and lib_SE_index > 0:
                                            count_duo_SE += 1
                                            g = seq_only[lib_E_index+6:lib_E_index+16]
                                            k = seq_only_SE[lib_SE_index+6:lib_SE_index+26]
                                            if g in libE_dict and k in libS_dict:
                                                count_duo_SE_both += 1
                                                perfect_oligo.write("%s^%s\t%s^%s\n" % (k,g,libS_dict[k],libE_dict[g]))
                                                fmatch.write("%s\t%s^%s\t%s^%s\tSE\n" % (record.name,k,g,libS_dict[k],libE_dict[g]))
                                            if g in libE_dict and k not in libS_dict:
                                                partial_oligo.write("%s\t%s\t%s\t-\t%f\t%f\n" % (k,g,libE_dict[g],float(lib_E_index), float(lib_SE_index)))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                                count_duo_SE_no_S += 1
                                                # matched_all_oligo.write("%s^%s\t%s^\n" % (g,k,libP_dict[g]))
                                            if g not in libE_dict and k in libS_dict:
                                                partial_oligo.write("%s\t%s\t-\t%s\t%f\t%f\n" % (k,g,libS_dict[k],float(lib_E_index), float(lib_SE_index)))
                                                SeqIO.write(record, partial_fastq, "fastq")
                                                count_duo_SE_no_E += 1
                                                # matched_all_oligo.write("%s^%s\t^%s\n" % (g,k,libA_dict[k]))
                                            if g not in libE_dict and k not in libS_dict:
                                                # matched_all_oligo.write("%s^%s\t^\n" % (g,k))
                                                no_oligo_match += 1
                                                count_duo_SE_no_oligo += 1
                                                SeqIO.write(record, no_match_oligo, "fastq")
                                        if check_S_index == -1 and check_E_index == 1:
                                            count_duo_no_match += 1
                                            count_duo_no_match_ES += 1
                                            count_duo_ES += 1
                                            SeqIO.write(record, no_match_duo, "fastq")
                                        if check_E_index == -1 and check_S_index == 1:
                                            count_duo_no_match += 1
                                            count_duo_no_match_SE += 1
                                            count_duo_SE += 1
                                            SeqIO.write(record, no_match_duo, "fastq")
                                        if check_S_index == -1 and check_E_index == -1:
                                            count_duo_no_match += -1
                                        # if lib_P_index == 0 and lib_A_index == 0 and lib_PiA_index == 0 and lib_AiP_index == 0:
                                        #     count_duo_no_match += 1
                                            SeqIO.write(record, no_match_duo, "fastq")
                                    if i%100 == 0:
                                        print(i)
                                    i += 1
                                    # if i > 550:
                                    #     break

print("%s seconds" % (time.time() - start_time))
with open ("%s_statistics.txt"  %  out_dir, "w") as run_stats:
    run_stats.write("Total Records Checked: %f \n" % float(i))
    run_stats.write("Percent Matched to Library: %f \n" % float(1-((count_duo_no_match+count_single_no_match)/i)))
    run_stats.write("Percent Duo Matched: %f \n" % float(1-(count_duo_no_match/count_duo)))
    run_stats.write("Percent Single Matched: %f \n" % float(1-(count_single_no_match/count_single)))
    run_stats.write("Percent Matched with No Oligos: %f\n" % float(no_oligo_match/((count_duo-count_duo_no_match)+(count_single-count_single_no_match))))
    run_stats.write("Percent S into E Mis-Match: %f\t%f\t%f\n" % (float(count_duo_no_match_ES/count_duo_ES), float(count_duo_no_match_ES), float(count_duo_ES)))
    run_stats.write("Percent E into S Mis-Match: %f\t%f\t%f\n" % (float(count_duo_no_match_SE/count_duo_SE), float(count_duo_no_match_SE), float(count_duo_SE)))
    run_stats.write("Single E: \n \tPercent Oligos Matched: %f\n \tPercent Oligos Not Matched: %f\n" % (float(1-(count_single_E_no_match/count_single_E)), float(count_single_E_no_match/count_single_E)))
    run_stats.write("Single S: \n \tPercent Oligos Matched: %f\n \tPercent Oligos Not Matched: %f\n" % (float(1-(count_single_S_no_match/count_single_S)), float(count_single_S_no_match/count_single_S)))
    run_stats.write("E into S: \n \tPercent Both Oligos Matched: %f\n \tPercent No Oligos Matched: %f\n \tPercent only E oligo matched: %f\n \tPercent only S oligo matched: %f\n" % (float(count_duo_SE_both/count_duo_SE), float(count_duo_SE_no_oligo/count_duo_SE), float(count_duo_SE_no_S/count_duo_SE), float(count_duo_SE_no_E/count_duo_SE)))
    run_stats.write("S into E: \n \tPercent Both Oligos Matched: %f\n \tPercent No Oligos Matched: %f\n \tPercent only E oligo matched: %f\n \tPercent only S oligo matched: %f\n" % (float(count_duo_ES_both/count_duo_ES), float(count_duo_ES_no_oligo/count_duo_ES), float(count_duo_ES_no_S/count_duo_ES), float(count_duo_ES_no_E/count_duo_ES)))
