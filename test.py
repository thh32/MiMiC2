import os

# Define output directory (default to current directory)
output_dir = os.getenv("OUTPUT_DIR", ".")

# Set the output file path
outfile = os.path.join(output_dir, "results", "testfile.txt")

# Write to the file
with open(outfile, 'w') as file:
    file.write("testing line1\n")
    file.write("testing line2\n")

