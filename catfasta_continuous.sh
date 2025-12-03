#!/bin/zsh

# remove any pre-existing fasta files
rm -f sample*.fasta

# for each sample out of A B C and D
for sample in A B C D; do

    # create new sample fasta file with header
    echo ">Sample${sample}" > sample${sample}_cont.fasta

    # create temporary working fasta to collect parts
    touch wrk${sample}.fasta
  
    # process each part (1, 2, 3) separately
    for part in 1 2 3; do
    
        # filter bases by phred score 20 (99% accuracy) using seqtk, substitute masked bases with N
        # remove spaces and sample headers using sed, extract sequence using grep
        # trim new line then append to working fasta file
        seqtk seq -q 20 sample${sample}*part${part}.FASTQ | sed 's/[atcgnyrwskmdvhbxn]/N/g' \
        | sed -e 's/ //g' -e 's/@sample._part_.//g' | grep '^[ATCGNYRWSKMDVHBXN-]*$' \
        | tr -d '\n' >> wrk${sample}.fasta

    done

    # use fold to limit line length to 80 characters (typical fasta format)
    fold -w 80 wrk${sample}.fasta >> sample${sample}_cont.fasta

done

# remove working files
rm -f wrk*.fasta