library(data.table)

# Path to the combined CSV file
combined_csv_file <- 'combined_p2x_filtered.csv'
output_csv_file <- 'combined_p2x_filtered_no_duplicates.csv'

# Read the header row to identify duplicate columns
header <- fread(combined_csv_file, nrows = 0, header = TRUE)

# Initialize a list to track original sample names and stripped sample names
original_names <- colnames(header)[-1]  # Exclude the first entry which is the row index name
stripped_names <- sub("_[^_]+_p2x_filtered$", "", original_names)

# Check for duplicates based on stripped sample names
duplicate_indices <- which(duplicated(stripped_names))

# Print some of the original and stripped names for verification
cat("Some of the original and stripped names:\n")
for (i in 1:min(10, length(original_names))) {
  cat(sprintf("Original: %s -> Stripped: %s\n", original_names[i], stripped_names[i]))
}

# Print duplicate columns with details
cat("Duplicate cells (column headers) found:\n")
for (i in duplicate_indices) {
  cat(sprintf("Duplicate: %s\n", original_names[i]))
}
cat(sprintf("Total number of duplicates: %d out of %d headers\n", length(duplicate_indices), length(original_names)))

# Read the entire CSV file
df <- fread(combined_csv_file, header = TRUE)

# Remove duplicate columns
df_no_duplicates <- df[, -c(duplicate_indices + 1), with = FALSE]  # +1 to account for the row index column

# Save the DataFrame without duplicates to a new CSV file
fwrite(df_no_duplicates, output_csv_file)

cat(sprintf("Duplicates removed. New CSV file saved as %s\n", output_csv_file))
