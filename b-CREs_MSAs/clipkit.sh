#!/bin/bash

# Find all folders ending with _Cleaned_MSAs in Intermediate_Files directory
INTERMEDIATE_DIR="Intermediate_Files"

# Check if Intermediate_Files directory exists
if [ ! -d "$INTERMEDIATE_DIR" ]; then
    echo "Error: $INTERMEDIATE_DIR directory not found"
    exit 1
fi

for INPUT_DIR in "$INTERMEDIATE_DIR"/*_Cleaned_MSAs; do
    # Skip if no matching folders found
    [ -d "$INPUT_DIR" ] || continue
    
    # Get the folder name without path
    FOLDER_NAME=$(basename "$INPUT_DIR")
    
    # Generate output directory name by replacing _Cleaned_MSAs with _MSAs_Ready
    OUTPUT_DIR_NAME="${FOLDER_NAME/_Cleaned_MSAs/_MSAs_Ready}"
    OUTPUT_DIR="$INTERMEDIATE_DIR/$OUTPUT_DIR_NAME"
    
    # Check if output directory already exists and has .fa files
    if [ -d "$OUTPUT_DIR" ] && [ -n "$(find "$OUTPUT_DIR" -maxdepth 1 -name "*.fa" -print -quit)" ]; then
        echo "Folder '$OUTPUT_DIR_NAME' already processed - skipping"
        continue
    fi
    
    # Number of parallel processes to run
    NUM_PROCESSES=7  # Adjust based on your CPU cores
    
    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR"
    
    # Function to process a single file
    process_file() {
        local file="$1"
        local output_dir="$2"
        
        # Extract the base name of the file (without extension)
        base_name=$(basename "$file" .fa)
        
        # Remove "cleaned_" prefix if present
        output_name=${base_name#cleaned_}
        
        # Run clipkit with kpic-gappy mode
        clipkit "$file" -m kpic-gappy -o "$output_dir/${output_name}.fa"
        
        # Print a message to indicate progress
        echo "Processed $file -> $output_dir/${output_name}.fa"
    }
    
    # Export the function so it's available to parallel
    export -f process_file
    
    echo ""
    echo "Processing $FOLDER_NAME..."
    
    # Find all .fa files and process them in parallel
    find "$INPUT_DIR" -name "*.fa" | xargs -P "$NUM_PROCESSES" -I{} bash -c "process_file '{}' '$OUTPUT_DIR'"
    
    echo "Completed processing $FOLDER_NAME -> $OUTPUT_DIR_NAME"
    echo ""
done

echo "All folders processed."