#!/bin/bash

# Script to change ".mzML" to ".raw" on line 7 of .mztab files in the current directory

# Directory containing the .mztab files
DIRECTORY="/storage/khangl/MS_benchmarking/Casanovo_output"

# Iterate over each .mztab file in the provided directory
for file in "$DIRECTORY"/*.mztab; do
    if [[ -f "$file" ]]; then
        echo "Processing file: $file"
        
        # Use 'sed' to perform the in-place replacement on line 7
        sed -i '7s/.mzML/.raw/' "$file"
        
        if [[ $? -eq 0 ]]; then
            echo "Successfully updated $file"
        else
            echo "Failed to update $file"
        fi
    else
        echo "No .mztab files found in the directory"
    fi
done

echo "Processing completed."