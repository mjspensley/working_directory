#!/bin/zsh

# remove any pre-existing fasta files
rm -f sample*.fasta

# for each sample out of A B C and D
for i in A B C D; do

  # create new sample fasta file with header
  echo ">sample${i}" > sample${i}.fasta

  # use sed to remove spaces and sample headers then grep to retrieve  sequence and append to working fasta file
  sed -e 's/ //g' -e 's/@sample._part_.//g' sample${i}*.FASTQ |  grep '^[ATCGNYRWSKMDVHBXN-]*$' >> wrk${i}.fasta

  # remove empty lines
  sed -i '' '/^$/d' wrk${i}.fasta

  # concatenate lines together to make one continuous sequence and append to final fasta file
  paste -sd '\0' wrk${i}.fasta >> sample${i}.fasta

done

# remove working fasta file
rm -f wrk*.fasta
