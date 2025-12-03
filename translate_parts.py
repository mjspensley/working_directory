# load glob
import glob

# load biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# match fasta files containing sequences to be translated
matched_fastas = glob.glob("sample?_parts.fasta")    # stores any fastas matched using the specified pattern

# iterate over matched fastas
for fasta in matched_fastas:

    # store sample and part identifiers
    sample = fasta[6]

    # open matched fasta as input file
    with open(fasta) as nt_sequences_fasta:

        # load sequences from input file
        nt_sequences = list(SeqIO.parse(nt_sequences_fasta, "fasta"))

    # create list to store output amino acid SeqRecords
    aa_sequences = []

    # loop over each nucleotide sequence record parsed from input fasta files
    for record in nt_sequences:

        # store the nucleotide sequence to be translated as a string
        dna_seq = record.seq

        # store the best open reading frame found and the frame in which it is contained
        best_orf = ""               # longest ORF found across all frames
        best_frame = None           # reading frame (0,1,2) containing the longest ORF

        # loop over nucleotide sequence to translate in all three forward reading frames
        for frame in range(3):

            # translate from this frame
            ## use vertebrate mitochondrial code for translation of codons (table=2)
            ## keep stop codons (to_stop=False)
            protein = dna_seq[frame:].translate(table=2, to_stop=False)
            
            # split translated sequence into fragments separated by stop codons ("*")
            fragments = str(protein).split("*")         # each fragment is a potential ORF
            
            # determine the longest amino acid sequence (fragment) in this frame
            frame_longest = max(fragments, key=len)
            
            # if the longest orf in the present frame is longer than any previously identified...
            if len(frame_longest) > len(best_orf):

                # overwrite and store in previously defined fields for best orf and reading frame
                best_orf = frame_longest
                best_frame = frame

        # convert longest orf string back into Seq object
        best_protein = Seq(best_orf)

        # parse the species name (Latin binomial) where possible from the FASTA description
        ## replace any spaces with underscores to prevent issues
        start = record.description.find("[organism=") + len("[organism=")
        end = record.description.find("]", start)
        latin_name = record.description[start:end].replace(" ", "_")

        # create a new SeqRecord containing the longest ORF for this sequence
        aa_record = SeqRecord(
            best_protein,
            id=record.id,                      # keeps the original FASTA ID
            name=latin_name,                   # stores the species name
            description=record.description,    # preserves the full original header
            annotations={                      # stores the reading frame from which the ORF was obtained
                "type": "longest ORF",
                "frame": best_frame            
            }
        )

        # append new SeqRecord to output list
        aa_sequences.append(aa_record)

    # write all longest ORF sequences to a new FASTA file
    SeqIO.write(aa_sequences, f"sample{sample}_translated.fasta", "fasta")