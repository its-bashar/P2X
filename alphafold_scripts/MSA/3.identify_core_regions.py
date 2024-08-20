from Bio import AlignIO
from Bio.Align import AlignInfo

# Step 1: Load the Clustal Omega output file
alignment_file = "clustalo-I20240610-231846-0948-74374821-p1m.aln-fasta"
alignment = AlignIO.read(alignment_file, "fasta")

# Step 2: Identify conserved regions
alignment_length = alignment.get_alignment_length()
num_sequences = len(alignment)

conserved_positions = []

for i in range(alignment_length):
    column = alignment[:, i]
    if column[0] != '-':  # Ignore gaps
        if all(residue == column[0] for residue in column):
            conserved_positions.append(i)

# Helper function to find the N and C termini
def find_terminals(conserved_positions):
    if conserved_positions:
        n_terminal = conserved_positions[0]
        c_terminal = conserved_positions[-1]
    else:
        n_terminal = None
        c_terminal = None
    return n_terminal, c_terminal

n_terminal, c_terminal = find_terminals(conserved_positions)

# Step 3: Print N and C Terminus Regions
for record in alignment:
    seq_id = record.id
    sequence = str(record.seq)

    if n_terminal is not None and c_terminal is not None:
        n_terminal_seq = sequence[:n_terminal]
        c_terminal_seq = sequence[c_terminal+1:]
    else:
        n_terminal_seq = ""
        c_terminal_seq = ""

    print(f"Receptor: {seq_id}")
    print(f"N-terminal (non-conserved): {n_terminal_seq}")
    print(f"C-terminal (non-conserved): {c_terminal_seq}")
    print()


