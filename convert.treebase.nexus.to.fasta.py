#!/usr/bin/env python3

# Convert a NEXUS file from TreeBASE into a FASTA file that lots of DNA analysis software can understand.

# Handles special characters well, at least the ones I have come across so far.   Older versions of this code didn't.

# By Alan Rockefeller - March 28, 2025


import sys
import re

def nexus_to_fasta(input_file, output_file):
    """
    Converts DNA sequences from a NEXUS file to FASTA format.
    
    This function reads a NEXUS file and extracts DNA sequences from its MATRIX section,
    using taxon names from the TAXLABELS section. Quoted taxon names are preserved and cleaned
    (by removing quotes and replacing spaces with underscores) to ensure compatibility with
    FASTA formatting. The extracted sequences are written to the output file in blocks of 60
    characters per line. Status messages are printed to indicate success or any issues encountered.
        
    Args:
        input_file (str): Path to the input NEXUS file.
        output_file (str): Path to the output FASTA file.
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Find the MATRIX section in the NEXUS file
        matrix_match = re.search(r'MATRIX\s*(.*?);', content, re.DOTALL)
        if not matrix_match:
            print("Error: Could not find the MATRIX section in the NEXUS file.")
            return
        
        matrix_content = matrix_match.group(1).strip()
        
        # Find all taxon names in the TAXLABELS section to ensure proper handling
        taxlabels_match = re.search(r'TAXLABELS\s*(.*?);', content, re.DOTALL)
        taxon_names = []
        if taxlabels_match:
            taxlabels_content = taxlabels_match.group(1).strip()
            # Extract each taxon name, preserving quoted names as a single entity
            raw_names = re.findall(r"'[^']*'|\S+", taxlabels_content)
            for name in raw_names:
                # Clean the name for FASTA format
                clean_name = name.replace("'", "").replace('"', "").replace(" ", "_")
                taxon_names.append((name, clean_name))
        
        # Extract sequences from the matrix
        # Extract sequences from the matrix
        sequences = []
        # Process matrix content line by line for more efficient extraction
        matrix_lines = matrix_content.strip().split('\n')
        taxon_to_sequence = {}

        # First pass: identify starting positions for each taxon
        for line in matrix_lines:
            line = line.strip()
            if not line:
                continue
            
            # Try to match the beginning of a line with a taxon name
            for taxon_original, _ in taxon_names:
                # Handle quoted and unquoted taxon names
                if (line.startswith(taxon_original) or 
                    line.startswith(f"'{taxon_original}'") or 
                    line.startswith(f'"{taxon_original}"')):
                    # Extract the sequence part (everything after the taxon name)
                    seq_part = line[len(taxon_original):].strip()
                    if taxon_original not in taxon_to_sequence:
                        taxon_to_sequence[taxon_original] = seq_part
                    else:
                        taxon_to_sequence[taxon_original] += seq_part
        
        for taxon_original, taxon_clean in taxon_names:
            if taxon_original in taxon_to_sequence:
                sequence = taxon_to_sequence[taxon_original]
                # Remove all whitespace
                sequence = re.sub(r'\s+', '', sequence)
                sequences.append((taxon_clean, sequence))
            else:
                print(f"Warning: No sequence found for taxon '{taxon_original}'")
        # Write sequences in FASTA format
        with open(output_file, 'w', encoding='utf-8') as f:
            for taxon_name, sequence in sequences:
                f.write(f">{taxon_name}\n")
                # Write sequence in blocks of 60 characters
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + "\n")
        
        print(f"Successfully converted {len(sequences)} sequences to FASTA format.")
    
    except Exception as e:
        print(f"Error: {e}")

def main():
    """
    Entry point for converting DNA sequences from NEXUS to FASTA format.
    
    Validates the command-line arguments and invokes the conversion process. If the 
    required input and output file paths are not provided, the script prints a usage 
    message and exits.
    """
    if len(sys.argv) != 3:
        print("Usage: python convert.treebase.nexus.to.fasta.py input_file.nexus output_file.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    nexus_to_fasta(input_file, output_file)

if __name__ == "__main__":
    main()
