#!/bin/bash

# Loop over the range 0 to 7
for i in {0..7}
do
    # Create the filename
    output_file="90Sr_${i}"
    input_file="output_${i}.root"
    # Run the command with the file
    ./convert "$input_file" "$output_file" 1

    # Check if the command was successful
    if [ $? -eq 0 ]; then
    echo "Successfully processed $output_file"
    else
    echo "Error processing $output_file"
    fi
done