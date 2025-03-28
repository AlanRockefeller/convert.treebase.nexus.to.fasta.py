#!/usr/bin/env python3

# Convert a NEXUS file from TreeBASE into a FASTA file that lots of DNA analysis software can understand.

# Handles special characters well, at least the ones I have come across so far.   Older versions of this code didn't.

# By Alan Rockefeller - March 28, 2025


import sys
import re

def nexus_to_fasta(input_file, output_file):
    """
    Converts DNA sequences from NEXUS format to FASTA format.
    
    Args:
        input_file (str): Path to the input NEXUS file
        output_file (str): Path to the output FASTA file
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
        sequences = []
        for taxon_original, taxon_clean in taxon_names:
            # Find the sequence for this taxon in the matrix content
            # Handle the case where the taxon name might be quoted in the matrix
            pattern = re.escape(taxon_original) + r'\s+(.*?)(?=\s+\S+\s+|\s*$)'
            seq_match = re.search(pattern, matrix_content, re.DOTALL)
            if seq_match:
                sequence = seq_match.group(1).strip()
                # Remove all whitespace from the sequence
                sequence = re.sub(r'\s+', '', sequence)
                sequences.append((taxon_clean, sequence))
        
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
    if len(sys.argv) != 3:
        print("Usage: python nexus_to_fasta.py input_file.nexus output_file.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    nexus_to_fasta(input_file, output_file)

if __name__ == "__main__":
    main()
