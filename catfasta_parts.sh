#!/bin/zsh

# remove any pre-existing fasta files
rm -f sample*.fasta

# for each sample out of A B C and D
for sample in A B C D; do

  # create new sample fasta file
  touch sample${sample}_parts.fasta
  
  # process each part (1, 2, 3) separately
  for part in 1 2 3; do

    # add header for this part
    echo ">sample${sample}_part${part}" >> sample${sample}_parts.fasta
    
    # filter bases by phred score 20 (99% accuracy) using seqtk, substitute masked bases with N
    # remove spaces and sample headers using sed, extract sequence using grep and append to sample fasta
    # use fold to limit line length to 80 characters (typical fasta format)
    seqtk seq -q 20 sample${sample}*part${part}.FASTQ | sed 's/[atcgnyrwskmdvhbxn]/N/g' \
    | sed -e 's/ //g' -e 's/@sample._part_.//g' | grep '^[ATCGNYRWSKMDVHBXN-]*$' \
    | fold -w 80 >> sample${sample}_parts.fasta

  done

  # remove empty lines
  sed -i '' '/^$/d' sample${sample}_parts.fasta

done