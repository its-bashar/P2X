#!/bin/bash

# Define the input and output files
input_file="expression-summary-full-03-11-24.csv"
output_file="filtered_p2x_genes.csv"

# List of Ensembl IDs for P2X genes
ensembl_ids=("ENSG00000108405" "ENSG00000187848" "ENSG00000109991" "ENSG00000135124" "ENSG00000083454" "ENSG00000099957" "ENSG00000089041")

# Join the Ensembl IDs into a single regex pattern
ensembl_pattern=$(IFS="|"; echo "${ensembl_ids[*]}")

# Initialize the output file with the header
head -n 1 "$input_file" > "$output_file"

# Process the file line-by-line and append matching rows to the output file
awk -v pattern="$ensembl_pattern" -F, '
NR > 1 && $0 ~ pattern {
    print $0 >> "'"$output_file"'"
}
' "$input_file"

echo "Filtering completed. The output file is saved as $output_file."

