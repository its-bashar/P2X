import csv
import re
import pandas as pd
import gc

# Path to the combined CSV file
combined_csv_file = 'combined_p2x_filtered.csv'

# Initialize a list to track original sample names and stripped sample names
original_names = []
stripped_names = []

# Dictionary to store the indices of duplicates and their originals
duplicate_indices = {}

# Read the header row
print(f"Loading header from file: {combined_csv_file}")
with open(combined_csv_file, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)

# Total number of headers (excluding the first entry which is the row index name)
total_headers = len(header) - 1

# Process each column header (skip the first entry because it's the row index name)
print("Processing header row for duplicates...")
for index, col in enumerate(header[1:], start=1):
    # Store the original name
    original_names.append(col)

    # Strip the last part of the sample name (everything after the last underscore before the dataset id)
    stripped_name = re.sub(r'_[^_]+_p2x_filtered$', '', col)
    stripped_names.append(stripped_name)

# Check for duplicates based on stripped sample names
seen_stripped_names = {}
duplicate_columns = []

for index, stripped in enumerate(stripped_names):
    if stripped in seen_stripped_names:
        original_index = seen_stripped_names[stripped]
        duplicate_indices[index] = original_index
    else:
        seen_stripped_names[stripped] = index

# Print some of the original and stripped names for verification
print("Some of the original and stripped names:")
for original, stripped in zip(original_names[:10], stripped_names[:10]):
    print(f"Original: {original} -> Stripped: {stripped}")

# Print duplicate columns with details
print("Duplicate cells (column headers) found:")
for duplicate_index, original_index in duplicate_indices.items():
    print(f"Duplicate: {original_names[duplicate_index]} from Original: {original_names[original_index]}")
print(f"Total number of duplicates: {len(duplicate_indices)} out of {total_headers} headers")
