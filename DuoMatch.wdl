workflow DUOmatch {
  File read_a
  File read_b
  File reference_fasta
  Int read_b_number
  Int read_len
  Int seq_min
  String working_directory
  String id_out
  String link_A_bc
  String link_P_bc
  String link_1_bc
  String link_2_bc
  String link_A_oligo
  String link_P_oligo
  String end_A_oligo
  String end_P_oligo
  String out_directory

  call Flash { input:
                  read_a=read_a,
                  read_b=read_b,
                  read_len=read_len,
                  id_out=id_out
                }
  call Pull_Barcodes { input:
                          #pull = pull,
                          working_directory=working_directory,
                          flashed=Flash.out,
                          read_b_number=read_b_number,
                          id_out=id_out,
                          seq_min=seq_min,
                          link_A_bc=link_A_bc,
                          link_P_bc=link_P_bc,
                          link_1_bc=link_1_bc,
                          link_2_bc=link_2_bc,
                          link_A_oligo=link_A_oligo,
                          link_P_oligo=link_P_oligo,
                          end_A_oligo=end_A_oligo,
                          end_P_oligo=end_P_oligo
                        }
  call Rearrange { input:
                        matched_barcodes=Pull_Barcodes.out1,
                        id_out=id_out
                      }
  call MiniMap { input:
                    reference_fasta=reference_fasta,
                    organized_fasta=Rearrange.out,
                    id_out=id_out
                  }
  call SAM2MPRA { input:
                      #sam=sam,
                      working_directory=working_directory,
                      sam_file=MiniMap.out1,
                      id_out=id_out
                    }
  call Sort { input:
                  MPRA_out=SAM2MPRA.out,
                  id_out=id_out
                }
  call Ct_Seq { input:
                    #count=count,
                    working_directory=working_directory,
                    sorted=Sort.out,
                    id_out=id_out
                  }
  call Parse { input:
                #  parse=parse,
                  working_directory=working_directory,
                  counted=Ct_Seq.out,
                  id_out=id_out
                }
  call preseq { input:
                 counted=Ct_Seq.out,
                 id_out=id_out
              }
  call qc_plot { input:
                    parsed=Parse.out,
                    working_directory=working_directory,
                    id_out=id_out
              }
  call relocate { input:
                    flashed=Flash.out,
                    matched=Pull_Barcodes.out1,
                    rejected=Pull_Barcodes.out2,
                    organized_fasta=Rearrange.out,
                    sam_file=MiniMap.out1,
                    map_log=MiniMap.out2,
                    MPRA_out=SAM2MPRA.out,
                    sorted=Sort.out,
                    counted=Ct_Seq.out,
                    parsed=Parse.out,
                    preseq_hist=preseq.hist,
                    preseq_res=preseq.res,
                    qc_plot=qc_plot.plots,
                    out_directory=out_directory
                  }
  }

task Flash {
  # Flashing raw fastq files together
  File read_a
  File read_b
  Int read_len
  String id_out
  command {
    flash2 -r ${read_len} -f 274 -s 20 -o ${id_out}.merged -t 10 ${read_a} ${read_b}
    }
  output {
    File out="${id_out}.merged.extendedFrags.fastq"
    }
  }
task Pull_Barcodes{
  File flashed
  Int read_b_number
  Int seq_min
  String working_directory
  String id_out
  String link_A_bc
  String link_P_bc
  String link_1_bc
  String link_2_bc
  String link_A_oligo
  String link_P_oligo
  String end_A_oligo
  String end_P_oligo

  command {
    perl ${working_directory}/map_barcodes_duo.pl ${flashed} ${read_b_number} ${id_out} ${link_A_bc} ${link_P_bc} ${link_1_bc} ${link_2_bc} ${link_A_oligo} ${link_P_oligo} ${end_A_oligo} ${end_P_oligo} ${seq_min}
    }
  output {
    File out1="${id_out}.match"
    File out2="${id_out}.reject"
    }
  }
task Rearrange {
  File matched_barcodes
  String id_out
  command <<<
    awk '{print ">"$1"&"$6"#"$3"\n"$4}' ${matched_barcodes} > ${id_out}.merged.match.enh.fa
    >>>
  output {
    File out="${id_out}.merged.match.enh.fa"
    }
  }
task MiniMap {
  File reference_fasta
  File organized_fasta
  String id_out
  command {
    minimap2 --for-only -Y --secondary=no -m 10 -n 1 -t 30 --end-bonus 12 -O 5 -E 1 -k 10 -2K50m --MD --eqx --cs=long -c -a ${reference_fasta} ${organized_fasta} > ${id_out}.merged.match.enh.sam 2> ${id_out}.merged.match.enh.log
    }
  output {
    File out1="${id_out}.merged.match.enh.sam"
    File out2="${id_out}.merged.match.enh.log"
    }
  }
task SAM2MPRA {
  #File sam
  File sam_file
  String working_directory
  String id_out
  command {
    perl ${working_directory}/SAM2MPRA.pl -C ${sam_file} ${id_out}.merged.match.enh.mapped
    }
  output {
    File out="${id_out}.merged.match.enh.mapped"
    }
  }
task Sort {
  File MPRA_out
  String id_out
  command {
    sort -S30G -k2 ${MPRA_out} > ${id_out}.merged.match.enh.mapped.barcode.sort
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.sort"
    }
  }
task Ct_Seq {
  #File count
  File sorted
  String working_directory
  String id_out
  command {
    perl ${working_directory}/Ct_seq.pl ${sorted} 2 4 > ${id_out}.merged.match.enh.mapped.barcode.ct
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct"
    }
  }
task Parse {
  #File parse
  File counted
  String working_directory
  String id_out
  command {
    perl ${working_directory}/parse_map.pl ${counted} > ${id_out}.merged.match.enh.mapped.barcode.ct.parsed
    }
  output {
    File out="${id_out}.merged.match.enh.mapped.barcode.ct.parsed"
    }
  }
task preseq {
 File counted
 String id_out
 command <<<
   awk '{ct[$4]++}END{for (i in ct)print i "\t" ct[i]}' ${counted} | sort -k1n > ${id_out}.merged.match.enh.mapped.barcode.ct.hist
   preseq lc_extrap -H ${id_out}.merged.match.enh.mapped.barcode.ct.hist -o ${id_out}.merged.match.enh.mapped.barcode.ct.hist.preseq -s 25000000 -n 1000 -e 1000000000
   >>>
 output {
   File hist="${id_out}.merged.match.enh.mapped.barcode.ct.hist"
   File res="${id_out}.merged.match.enh.mapped.barcode.ct.hist.preseq"
   }
 }
task qc_plot {
  File parsed
  String working_directory
  String id_out
  command {
    Rscript ${working_directory}/mapping_QC_plots.R ${parsed} ${id_out}
    }
  output {
    File plots="${id_out}_barcode_qc.pdf"
    }
}
task relocate{
  File flashed
  File matched
  File rejected
  File organized_fasta
  File sam_file
  File map_log
  File MPRA_out
  File sorted
  File counted
  File parsed
  File preseq_hist
  File preseq_res
  File qc_plot
  String out_directory
  command {
      mv ${flashed} ${matched} ${rejected} ${organized_fasta} ${sam_file} ${map_log} ${MPRA_out} ${sorted} ${counted} ${parsed} ${preseq_hist} ${preseq_res} ${qc_plot} ${out_directory}
    }
  }
