#!/bin/bash

# Define input and output directories here
INPUT_DIR="Intermediate_Files/cleaned"
OUTPUT_DIR="Neutral"

# Number of parallel processes to run
NUM_PROCESSES=32  # Adjust based on your CPU cores

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to process a single file
process_file() {
    local file="$1"
    
    # Extract the base name of the file (without extension)
    base_name=$(basename "$file" .fa)
    
    # Remove "cleaned_" prefix if present
    output_name=${base_name#cleaned_}
    
    # Run clipkit with kpi-gappy mode and gap threshold of 0.2
    clipkit "$file" -m kpic-gappy -g 0.2 -o "$OUTPUT_DIR/${output_name}.fa"
    
    # Print a message to indicate progress
    echo "Processed $file -> $OUTPUT_DIR/${output_name}.fa"
}

# Export the function and directories so they're available to parallel
export -f process_file
export OUTPUT_DIR

# Find all .fa files and process them in parallel
find "$INPUT_DIR" -name "*.fa" | xargs -P "$NUM_PROCESSES" -I{} bash -c 'process_file "{}"'

echo "All files processed. Results saved in $OUTPUT_DIR"