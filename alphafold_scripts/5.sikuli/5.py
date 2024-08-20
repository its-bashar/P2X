# Import necessary libraries
import os
import itertools

# Define the directory path to sequence files
sequence_dir = os.path.join(os.path.dirname(getBundlePath()), "Sequences")

# Define the sequence files for all P2X receptors
sequence_files = {
    "P2X1": "P2RX1_HUMAN.txt",
    "P2X2": "P2RX2_HUMAN.txt",
    "P2X3": "P2RX3_HUMAN.txt",
    "P2X4": "P2RX4_HUMAN.txt",
    "P2X5": "P2RX5_HUMAN.txt",
    "P2X6": "P2RX6_HUMAN.txt",
    "P2X7": "P2RX7_HUMAN.txt"
}

# Function to read the sequence from a file
def read_sequence(file_path):
    with open(file_path, 'r') as file:
        return file.read().strip()

# Function to automate the sequence submission process
def submit_sequence(seq1_name, seq2_name, seq3_name, seq1, seq2, seq3):
    # Clear previous inputs
    click("1718216377190.png")  # Embed the screenshot of the Clear button
    wait(1)
    # Paste the first sequence
    click(Pattern("1718216526701.png").targetOffset(-25,-25))  # Embed the screenshot of the Paste Sequence box
    paste(seq1)
    
    # Paste the second sequence
    click("1718228184249.png") # Embed the screenshot of "add entity" with the space next to it (make it press on the space next to it)
    wait(1)
    click("1718216573086.png")  # Embed the screenshot of the Add Entity button
    click(Pattern("1718216526701.png").targetOffset(-25,-25))  # Embed the screenshot of the Paste Sequence box
    paste(seq2)
    
    # Paste the third sequence
    click("1718228184249.png")# Embed the screenshot of "add entity" with the space next to it (make it press on the space next to it)

    wait(1)
    click("1718216573086.png")  # Embed the screenshot of the Add Entity button
    click(Pattern("1718216526701.png").targetOffset(-25,-25))  # Embed the screenshot of the Paste Sequence box
    paste(seq3)

    # To Test uncomment this and comment the uncommented line below it
    #click("1718216594654.png")  # Embed the screenshot of the Continue and Preview Job button
    #click("1718225301287.png")

    click("1718228184249.png")# Embed the screenshot of "add entity" with the space next to it (make it press on the space next to it)

    wait(1)
    
    # add oliec acid OLA
    #click("1718216573086.png")  # Embed the screenshot of the Add Entity button
    #click(Pattern("1718216526701.png").targetOffset(-372,-20))  # Embed the molecule type
    #click(Pattern("1720816346853.png").targetOffset(-7,72))  # Embed "Ligand" choice
    #click("1720816474623.png")  # Embed copies choice
    #type("a", KeyModifier.CTRL)
    #paste("50")
    #click("1720816548628.png")  # Embed ligand choices
    
    # Press the down button 12 times
    #for i in range(12):
    #    type(Key.DOWN)
    #    wait(0.5)  # Wait for 0.5 seconds between each press
    # Press Enter
    #type(Key.ENTER)
    #click("1718228184249.png")# Embed the screenshot of "add entity" with the space next to it (make it press on the space next to it)
    #wait(1)
    click("1718226675606.png")# Embed the screenshot of "Continue and preivew job"

    click("1718226487045.png")# Embed the screenshot of part of the file name (do NOT screenshot the deualt file name as it changes for each submmition)
    type("a", KeyModifier.CTRL)
    paste("{}_{}_{}_TRIMMED_SPLICE_VARIENT".format(seq1_name, seq2_name, seq3_name))  # Write the file name

    # setting the seed
    click(Pattern("1718227129754.png").targetOffset(13,0))# Embed the screenshot of seed toogle
    click("1718227178314.png")# Embed the screenshot of seed number
    type("a", KeyModifier.CTRL)
    paste("1123581321")
    click("1718227274353.png")# Embed the screenshot of "Confirm and submit job"
    
# Generate all unique combinations for all P2X receptors (considering order doesn't matter)
p2x_receptors = ['P2X1', 'P2X2', 'P2X3', 'P2X4', 'P2X5', 'P2X6', 'P2X7']
combinations = set(itertools.combinations_with_replacement(p2x_receptors, 3))

# Convert combinations to the required format and sort them
formatted_combinations = ['-'.join(combination) for combination in combinations]
formatted_combinations.sort()

# Function to get sequences to be submitted from files
def get_sequences(seq_combination):
    return [seq_combination, read_sequence(os.path.join(sequence_dir, sequence_files[seq_combination[0]])), 
            read_sequence(os.path.join(sequence_dir, sequence_files[seq_combination[1]])), 
            read_sequence(os.path.join(sequence_dir, sequence_files[seq_combination[2]]))]

# Log file path
log_file_path = os.path.join(os.path.dirname(getBundlePath()), "submission_log.txt")

# Function to read the log file
def read_log():
    if not os.path.exists(log_file_path):
        return []
    with open(log_file_path, 'r') as file:
        return file.read().splitlines()

# Function to write to the log file
def write_log(combination):
    with open(log_file_path, 'a') as file:
        file.write(combination + '\n')

# Read completed combinations from log file
completed_combinations = read_log()

# Filter out completed combinations from the sequences list
remaining_combinations = [seq for seq in formatted_combinations if seq not in completed_combinations]

# Limit to 20 submissions
combinations_to_submit = remaining_combinations[:20]

# Test the script with the filtered sequences
for combination in combinations_to_submit:
    seq_combination = combination.split('-')
    sequences = get_sequences(seq_combination)
    submit_sequence(seq_combination[0], seq_combination[1], seq_combination[2], sequences[1], sequences[2], sequences[3])
    write_log(combination)
