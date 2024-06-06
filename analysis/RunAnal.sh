#!/bin/bash

# Define the file path
FILE="../build/output_files/output_t0.root"

# Check if the file exists
if [ -f "$FILE" ]; then
    echo "File '$FILE' exists. Moving to current directory."
    mv "$FILE" .
    if [ $? -eq 0 ]; then
        echo "File successfully moved to the current directory."
    else
        echo "Failed to move the file."
    fi
else
    echo "File '$FILE' does not exist."
fi

root 'RecoTrack.C("output_t0.root")'
echo "Executing Study..."
./study elab_output_t0.root -l