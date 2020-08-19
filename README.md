# MPRAduo

*An application of Massively Parallel Reporter Assays (MPRA) to determine the effect of silencers within the genome.*

## Before running the pipeline
* Have the latest version of Cromwell and Womtool in your workspace
  * `conda install -c bioconda cromwell`
  * `conda install -c bioconda womtool`

* Have modules for FLASH2, minimap2, preseq, pandas, and Biopython available
  * `conda install -c bioconda flash2 `
  * `conda install -c bioconda minimap2`
  * `conda install -c bioconda preseq`
  * `conda install -c anaconda pandas`
  * `conda install -c anaconda biopython`

* Make sure all the available scripts (except for the WDL itself) are in a known directory (you will need to provide the path to this directory)

## Running the WDL
* Validate the file
  `womtool validate <pipeline_name>.wdl`

  **NB: use the version number for your version of Womtool downloaded above**

* Generate inputs file
  `womtool inputs <pipeline_name>.wdl > <your_projects_name>_inputs.json`

  **NB: see the "Filling in the json" section below for detailed description of input needed**

* Run the pipeline with your inputs
  `cromwell run <pipeline_name>.wdl --inputs <your_projects_name>_inputs.json`

The set of pipelines and scripts outlined below are used for matching the oligos and barcodes, counting the sequenced duo barcodes, and a set of functions in R that can be used to effectively analyze the count output.

## WDL Pipeline Descriptions

### _DuoMatch.wdl_

![Grpaphical Pipeline](graphics/DuoMatch_pipeline.png?raw=true "DuoMatch Graphical Pipeline")

Takes read 1 and read 2 of the barcode-oligo sequences, flashes them together (`FLASH2`), and pulls the barcode, oligo and UMI (`map_barcodes_duo.pl`). These sequences are then reorganized into a fasta format, with the barcode and UMI tacked on to the record ID and the oligo sequence as the sequence, and then mapped (`minimap2`) using the oligo order fasta as the reference. The matching oligo ID, barcode, CIGAR information, among other data points are pulled from the SAM output (`SAM2MPRA.pl`), and then counted (`Ct_seq.pl`) and parsed (`parse_map.pl`). `preseq` is also run, for the purposes of use in QC.

An example inputs file can be found in this repository, and a description of each of the inputs required can be found below:

```
    {
      "DuoMatch.read_a": "/full/path/to/read/1.fastq.gz",
      "DuoMatch.read_b": "/full/path/to/read/2.fastq.gz",
      "DuoMatch.reference_fasta": "/full/path/to/oligo/order/fasta.fa",
      "DuoMatch.read_b_number": "2 (unless read b above is read 1)",
      "DuoMatch.seq_min": "100 (suggested, minimum length of barcode-oligo sequence)",
      "DuoMatch.working_directory": "/path/to/directory/of/scripts",
      "DuoMatch.id_out": "Dictionary_ID_Name",
      "DuoMatch.link_A_bc": "TCTAGA (6 bases at the barcode end of the linker sequence for Silencers)",
      "DuoMatch.link_P_bc": "AGTCAG (6 bases at the barcode end of the linker sequence for Enhancers)",
      "DuoMatch.link_1_bc": "CTAG (4 bases immediately following the UMI - Enhancers)",
      "DuoMatch.link_2_bc": "CTCA (4 bases immediately following the UMI - Silencers)",
      "DuoMatch.link_A_oligo": "AGTG (4 bases immediately preceding the oligo - Silencers)",
      "DuoMatch.link_P_oligo": "AGGG (4 bases immediately preceding the oligo - Enhancers)",
      "DuoMatch.end_A_oligo": "CGTC (4 bases immediately following the oligo - Silencers)",
      "DuoMatch.end_P_oligo": "GCAA (4 bases immediately following the oligo - Enhancers)"
    }
```

**NB** The orientation for the sequences in order to find the linker indicators for the input file the UMI end should be on the left and the end of the oligo should be on the right.

#### Files to be used later:

Once the `DuoMatch` pipeline is run, the final output file will need to be passed to the `DuoCount` pipeline. Look in the first ~30 lines of the output to find the cromwell id for the relevant run. The file name ends with `.parsed` and can be found within the following path: `cromwell_executions/DuoMatch/<cromwell_id>/call-Parse/execution/<id_out>.merged.match.enh.mapped.barcode.ct.parsed`. The pipeline should be run twice, once for silencers and once for enhancers, and both parsed files should be passed to the `DuoCount` pipeline.

### _DuoCount.wdl_ and _DuoCountSingle.wdl_

![Graphical Pipeline](graphics/DuoCount_pipeline.png?raw-true "DuoCount Graphical Pipeline")

These two pipelines are virtually identical, the only difference is that `DuoCount.wdl` is for paired end reads and `DuoCountSingle.wdl` is for single end reads, both are expecting ~150bp reads. Read 1 (and read 2, if paired end reads) fastq files and the associated replicate IDs are passed to the pipeline (`library_separation_FLASH2_input.py`) along with the dictionary files created by `DuoMatch.wdl` to have the libraries sorted; single and duo libraries are separated and then identified as E,S,ES, or SE. The matched file, containing the record name, barcode, matched oligo, and the assigned library, for each replicate is formatted (`associate_tags_duo.pl`) and passed to a script (`compile_bc_duo.pl`) to assemble the information into a barcode level count table.

An example of the inputs file for each version of the script can be found in the repository, an example of the paired end inputs file can be found below:

```
    {
      "ReplicateCount.parsed_E": "/path/to/enhancers/dictionary/file",
      "ReplicateCount.parsed_S": "/path/to/silencers/dictionary/file",
      "ReplicateCount.read_a": ["/path/to/R1/celltype1/rep1.fastq.gz","/path/to/R1/celltype1/rep2.fastq.gz","/path/to/R1/celltype1/rep3.fastq.gz",...,"/path/to/R1/celltypen/repx.fastq.gz"],
      "ReplicateCount.read_b": ["/path/to/R2/celltype1/rep1.fastq.gz","/path/to/R2/celltype1/rep2.fastq.gz","/path/to/R2/celltype1/rep3.fastq.gz",...,"/path/to/R2/celltypen/repx.fastq.gz"], (NB: only required for paired end reads)
      "ReplicateCount.replicate_id": ["Celltype1_r1","Celltype1_r2","Celltype1_r3",..."Celltypen_rx"],
      "ReplicateCount.read_b_number": "2 (unless for paired read_b above is R1)",
      "ReplicateCount.id_out": "Overall_Project_ID",
      "ReplicateCount.working_directory": "/path/to/directory/of/scripts"
    }
```

#### Files to be used later:

For the following files the barcode-oligo format is as follows for duo libraries:
  * ES: enhancer^silencer
  * SE: silencer^enhancer

There are several files that will be useful in further analysis.
  * `library_separation_FLASH2_input.py` outputs found in `cromwell-executions/ReplicateCount/<cromwell_id/call-prep_counts/shard-*/execution/`:
    * `<replicate_id>.match`                  : tab separated file containing one record per line. The columns represent sequence ID, barcode, oligo, and library assignment.
    * `<replicate_id>_single_match.fasta`     : multiline fasta of all records identified as belonging to a single library.
    * `<replicate_id>_perfect_match.txt`      : tab separated file of perfect duo matches (both barcodes are present in the dictionary).
    * `<replicate_id>_partial_match.fastq`    : fastq file of all the records which had only one barcode present in the dictionary.
    * `<replicate_id>_partial_match.txt`      : tab separated file with columns for each barcode, oligo and the location the linker was found.
    * `<replicate_id>_single_no_match.fastq`  : fastq file of all the records, identified as single, where there was no linker sequence present.
    * `<replicate_id>_duo_no_match.fastq`     : fastq file of all the records, identified as duo, where there were no linker sequences present.
    * `<replicate_id>_oligo_no_match.fastq`   : fastq file of all the records, where linkers were found but barcodes were not in the dictionary.
    * `<replicate_id>_statistics.txt` :
      * Total Records checked
      * Total Records GFP not found
      * Percent matched to library
      * Percent Duo matched
      * Percent Single matched
      * Percent matched with no oligos
      * Percent ES Mis-Match
      * Percent SE Mis-Match
      * Single E:
        * Percent Oligos Matched
        * Percent Oligos not Matched
      * Single S:
        * Percent Oligos Matched
        * Percent Oligos not Matched
      * SE:
        * Percent Both Oligos Matched
        * Percent No Oligos Matched
        * Percent only E oligo Matched
        * Percent only S Oligo Matched
      * ES:
        * Percent Both Oligos Matched
        * Percent No Oligos Matched
        * Percent only E oligo Matched
        * Percent only S Oligo Matched
  * `<id_out>.count` found in `cromwell-executions/ReplicateCount/<cromwell_id>/call-make_count_table/execution/`:
    * Tab separated count file with columns: Barcode, Oligo, library, Replicate Ids
    * This file will then be used with the R functions described below.
