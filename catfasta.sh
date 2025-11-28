#!/bin/zsh

# remove any pre-existing fasta files
rm -f sample*.fasta

# for each sample out of A B C and D
for sample in A B C D; do

  # create new sample fasta file
  echo ">sample${sample}" >> sample${sample}.fasta
  
  # filter bases by phred score 20 (99% accuracy) using seqtk, substitute masked bases with N
  # remove spaces and sample headers using sed, extract sequence using grep and append to working fasta
  seqtk seq -q 20 sample${sample}*.FASTQ | sed 's/[atcgnyrwskmdvhbxn]/N/g' | sed -e 's/ //g' -e 's/@sample._part_.//g' | grep '^[ATCGNYRWSKMDVHBXN-]*$' >> wrk${sample}.fasta

  # remove empty lines from working fasta file
  sed -i '' '/^$/d' wrk${sample}.fasta

  # concatenate lines together to make one continuous sequence
  # create new line after 60 characters (fasta format) and append to final fasta file
  paste -sd '\0' wrk${sample}.fasta | fold -w 60 wrk${sample}.fasta >> sample${sample}.fasta

done

# remove working fasta
rm -f wrk*.fasta