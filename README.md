# MPRAduo

library_separation_FLASH2_input.py:
      Requires input of a FLASHed fastq file. Produces several files based on the records of the fastq.
      perfect_match.txt     - outputs list of barcodes and oligos that match perfectly based on the dictionary being used
      partial_match.fastq   - outputs the full record in fastq format of sequences which match to the library, but only find one oligo in the dictionary
      partial_match.txt     - outputs list of barcodes and oligos that partially match records based on the dictionary being used
      single_no_match.fastq - outputs the full record in fastq format of single library match where the linker can't be found
      duo_no_match.fastq    - outputs the full record in fastq format of duo library matches where no linker sequence can be found
      oligo_no_match.fastq  - outputs the full record in fastq format of sequences where linkers were matched but oligos couldn't be identified based on the dictionary being used.
      statistics.txt        - outputs various statistics used for quality control of library separation
