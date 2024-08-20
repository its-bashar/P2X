import csv
import os

input_file = 'combined_p2x_filtered_no_duplicates.csv'
output_file = 'final_count_matrix.csv'

# Create the output file and open it for writing
with open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)

    # Open the input file and process it line by line
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile)
        
        # Read all rows into a list
        all_rows = list(reader)
        
        # Get the number of columns from the first row
        num_cols = len(all_rows[0])
        
        # Transpose the rows and columns
        for col_index in range(num_cols):
            transposed_row = [all_rows[row_index][col_index] for row_index in range(len(all_rows))]
            writer.writerow(transposed_row)

