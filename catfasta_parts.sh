#!/bin/zsh

# remove any pre-existing fasta files
rm -f sample*.fasta wrk*.fasta

# for each sample out of A B C and D
for i in A B C D; do

  # create new sample fasta file
  touch sample${i}.fasta

  # process each part (1, 2, 3) separately
  for part in 1 2 3; do

    # add a header for each part
    echo ">sample${i}_part${part}" >> sample${i}.fasta

    # use sed to remove spaces and sample headers then grep to retrieve sequence and append to working fasta file
    sed -e 's/ //g' -e 's/@sample._part_.//g' sample${i}*part${part}.FASTQ | grep '^[ATCGNYRWSKMDVHBXN-]*$' >> wrk${i}_part${part}.fasta

    # remove empty lines
    sed -i '' '/^$/d' wrk${i}_part${part}.fasta

    # append working fasta file contents to sample fasta file
    cat wrk${i}_part${part}.fasta >> sample${i}.fasta

  done

done

# remove working fasta file
rm -f wrk*.fasta
