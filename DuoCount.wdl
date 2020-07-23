# Pipeline for counting MPRA Replicates
# Requires the same read_b_number as MPRAMatch
# Requires the parsed file output from MPRAMatch

workflow ReplicateCount {
  Array[File] read_a
  Array[File] read_b
  Array[String] replicate_id
  Array[Pair[File,File]] to_flash = zip(read_a,read_b)
  Array[Pair[String,Pair[File,File]]] id_flash = zip(replicate_id,to_flash)
  File parsed_S
  File parsed_E
  Int read_b_number
  String flags
  String id_out
  String working_directory

  scatter (replicate in id_flash) {
    Pair[File,File] reads = replicate.right
    call Flash { input:
                    read_a=reads.left,
                    read_b=reads.right,
                    id_out=replicate.left
                  }
    call prep_counts { input:
                          working_directory=working_directory,
                          sample_fastq=Flash.out,
                          parsed_S=parsed_S,
                          parsed_E=parsed_E,
                          read_b_number=read_b_number,
                          sample_id=replicate.left
                        }
    call associate { input:
                        working_directory=working_directory,
                        matched=prep_counts.out,
                        parsed=parsed,
                        read_b_number=read_b_number,
                        sample_id=replicate.left
                      }
                    }
  call make_infile { input:
                        working_directory=working_directory,
                        tag_files=associate.outF,
                        tag_ids=associate.outS,
                        id_out=id_out
                      }
  call make_count_table { input:
                            working_directory=working_directory,
                            list_inFile=make_infile.out,
                            flags=flags,
                            id_out=id_out
                          }
}

task Flash {
  # Flashing raw fastq files together
  File read_a
  File read_b
  String id_out
  command {
    flash2 -r 150 -f 274 -s 20 -o ${id_out}.merged -t 10 ${read_a} ${read_b}
    }
  output {
    File out="${id_out}.merged.extendedFrags.fastq"
    }
  }

task prep_counts {
  #File make_counts
  File sample_fastq
  File parsed_E
  File parsed_S
  Int read_b_number
  String working_directory
  String sample_id
  command {
    python ${working_directory}/library_separation_FLASH2_input.v2.py ${sample_fastq} ${parsed_E} ${parsed_S} ${sample_id}
    }
  output {
    File out="${sample_id}.match"
    }
  }
task associate {
  #File assoc
  File matched
  File parsed
  Int read_b_number
  String working_directory
  String sample_id
  command {
    perl ${working_directory}/associate_tags.pl ${matched} ${parsed} ${sample_id}.tag ${read_b_number}
    }
  output {
    File outF="${sample_id}.tag"
    String outS="${sample_id}"
    }
  }
task make_infile {
  Array[File] tag_files
  Array[String] tag_ids
  String working_directory
  String id_out
  command {
    python ${working_directory}/make_infile.py ${sep=',' tag_ids} ${sep=',' tag_files} ${id_out}
  }
  output {
    File out="${id_out}_samples.txt"
    }
  }
task make_count_table {
  #File compile
  File list_inFile
  String working_directory
  String? flags = ""
  String id_out
  command {
    perl ${working_directory}/compile_bc.pl ${flags} ${list_inFile} ${id_out}.count
    }
  output {
    File out="${id_out}.count"
    }
  }
