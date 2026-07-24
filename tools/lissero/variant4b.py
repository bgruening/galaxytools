import sys

def identify_variants_with_genes(input_file_path, output_file_path, simple_name):
    # Define the genes of interest
    genes_of_interest = ['LMO0737', 'ORF2110', 'ORF2819']

    # Open the input file and read its lines
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    # Check if the file has more than just the header
    if len(lines) <= 1:
        print("Input file does not contain enough data.")
        return

    # Extract the column headers and find the indices of the genes of interest
    headers = lines[0].strip().split('\t')
    gene_indices = [headers.index(gene) for gene in genes_of_interest]
    serotype_index = headers.index('SEROTYPE')

    # Modify the header to include the new first column
    modified_header = f"FileName\t{lines[0]}"
    
    # Initialize a list to hold the modified lines, starting with the modified header
    modified_lines = [modified_header]

    # Process each data line in the input file
    for line in lines[1:]:
        data = line.strip().split('\t')
        # Check if the genes of interest are all present (marked as "FULL")
        if all(data[index] == 'FULL' for index in gene_indices):
            # Modify the SEROTYPE column to "4b variant"
            data[serotype_index] = "4b variant"
        # Prepend the simple name to the line
        modified_line = f"{simple_name}\t" + '\t'.join(data) + '\n'
        # Add the modified line to the list
        modified_lines.append(modified_line)

    # Write the modified lines to the output file
    with open(output_file_path, 'w') as file:
        file.writelines(modified_lines)

    print(f'Results written to {output_file_path}')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file_path> <output_file_path> <simple_name>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    simple_name = sys.argv[3]
    identify_variants_with_genes(input_file_path, output_file_path, simple_name)

