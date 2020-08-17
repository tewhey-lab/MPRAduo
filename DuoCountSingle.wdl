# Pipeline for counting MPRA Replicates
# Specifically for singled end reads of replicates
# Requires the same read_b_number as MPRAMatch
# Requires the parsed file output from MPRAMatch

workflow ReplicateCount {
  Array[File] read_a
  Array[String] replicate_id
  Array[Pair[String,File]] id_rep = zip(replicate_id,read_a)
  File parsed_S
  File parsed_E
  Int read_b_number
  String id_out
  String working_directory

  scatter (replicate in id_rep) {
    call prep_counts { input:
                          working_directory=working_directory,
                          sample_fastq=replicate.right,
                          parsed_S=parsed_S,
                          parsed_E=parsed_E,
                          sample_id=replicate.left
                        }
    call associate { input:
                        working_directory=working_directory,
                        matched=prep_counts.out,
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
                            id_out=id_out
                          }
}

task prep_counts {
  #File make_counts
  File sample_fastq
  File parsed_E
  File parsed_S
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
  Int read_b_number
  String working_directory
  String sample_id
  command {
    perl ${working_directory}/associate_tags_duo.pl ${matched} ${sample_id}.tag ${read_b_number}
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
  String id_out
  command {
    perl ${working_directory}/compile_bc_duo.pl ${list_inFile} ${id_out}.count
    }
  output {
    File out="${id_out}.count"
    }
  }
