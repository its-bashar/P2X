import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Input file name
input_file = "P2X_receptors_trimmed.fasta"

# Directory to save the cleaned sequences
output_dir = "Sequences"

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Process the input file and clean the sequences
cleaned_sequences = []

with open(input_file, "r") as infile:
    # Parse the input file using SeqIO
    for record in SeqIO.parse(infile, "fasta"):
        # Remove hyphens from the sequence
        cleaned_seq = str(record.seq).replace("-", "")
        
        # Create a new SeqRecord with the cleaned sequence
        cleaned_record = SeqRecord(Seq(cleaned_seq), id=record.id, description=record.description)
        cleaned_sequences.append(cleaned_record)

# Save each cleaned sequence to a separate text file in the Sequences directory
for record in cleaned_sequences:
    # Get the sequence ID in the desired format (e.g., P2RX2_HUMAN)
    seq_id = record.id.split('|')[2].split('/')[0]  # Extract the part after the second '|' and before the first '/'
    
    # Create the file name
    file_name = os.path.join(output_dir, f"{seq_id}.txt")
    
    # Write the sequence to a text file
    with open(file_name, "w") as outfile:
        outfile.write(str(record.seq))

print("Cleaned sequences have been saved to individual text files in the 'Sequences' folder.")

