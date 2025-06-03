#!/usr/bin/env python3

# Convert a NEXUS file from TreeBASE into a FASTA file that lots of DNA analysis software can understand.

# Handles special characters well, at least the ones I have come across so far.   Older versions of this code didn't.

# By Alan Rockefeller - June 3, 2025

# Version 1.0


import sys
import re
import os

# Compiled regex patterns for use in nexus_to_fasta
MATRIX_PATTERN = re.compile(r'MATRIX\s*(.*?);', re.DOTALL)
TAXLABELS_PATTERN = re.compile(r'TAXLABELS\s*(.*?);', re.DOTALL)
TAXON_NAME_SPLIT_PATTERN = re.compile(r"'[^']*'|\S+")
SEQUENCE_WHITESPACE_PATTERN = re.compile(r'\s+')

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
    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
        return # Add return to stop execution in mocked environment
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
        return # Add return here too

    # Find the MATRIX section in the NEXUS file
    matrix_match = MATRIX_PATTERN.search(content)
    if not matrix_match:
        print("Error: Could not find the MATRIX section in the NEXUS file.")
        sys.exit(1)
        return # Add return to stop execution in mocked environment
    
    matrix_content = matrix_match.group(1).strip()
    
    try:
        # Find all taxon names in the TAXLABELS section to ensure proper handling
        taxlabels_match = TAXLABELS_PATTERN.search(content)
        taxon_names = []
        if taxlabels_match:
            taxlabels_content = taxlabels_match.group(1).strip()
            # Extract each taxon name, preserving quoted names as a single entity
            raw_names = TAXON_NAME_SPLIT_PATTERN.findall(taxlabels_content)
            for name_in_label in raw_names:
                # Store the original name (potentially quoted) for cleaning for FASTA header
                # Store the unquoted version as 'taxon_original' for matching in MATRIX
                
                original_name_for_matching = name_in_label
                if (name_in_label.startswith("'") and name_in_label.endswith("'")) or \
                   (name_in_label.startswith('"') and name_in_label.endswith('"')):
                    if len(name_in_label) > 1: # Avoid issues with just "'" or '"'
                        original_name_for_matching = name_in_label[1:-1]
                
                # Revert to simpler cleaning for FASTA header to isolate the issue
                # This was the logic that seemed most likely to be correct.
                clean_name_for_fasta = name_in_label.replace("'", "").replace('"', "").replace(" ", "_")
                
                taxon_names.append((original_name_for_matching, clean_name_for_fasta))
        
        # Refactored sequence extraction logic
        sequences = []
        matrix_lines = matrix_content.strip().split('\n')
        taxon_to_sequence_parts = {} # 1. Initialize dictionary

        # 2. Single pass through matrix_lines
        for line in matrix_lines:
            line = line.strip()
            if not line:
                continue
            
            # 3a. Determine which taxon this line belongs to
            # found_taxon_for_line = False # Keep for debugging if needed
            for taxon_original, _ in taxon_names:
                # 3b. Check for original_name (and quoted versions)
                prefix_to_check = ""
                name_len_in_line = 0
                actual_taxon_name_in_line = ""


                # Check quoted versions first
                quoted_single = f"'{taxon_original}'"
                quoted_double = f'"{taxon_original}"'

                if line.startswith(quoted_single):
                    prefix_to_check = quoted_single
                    name_len_in_line = len(quoted_single)
                    actual_taxon_name_in_line = taxon_original
                elif line.startswith(quoted_double):
                    prefix_to_check = quoted_double
                    name_len_in_line = len(quoted_double)
                    actual_taxon_name_in_line = taxon_original
                elif line.startswith(taxon_original):
                    # Check unquoted version last to avoid issues if a taxon name is a prefix of another's quoted version
                    prefix_to_check = taxon_original
                    name_len_in_line = len(taxon_original)
                    actual_taxon_name_in_line = taxon_original
                
                if prefix_to_check: # If any match type was found
                    # Ensure it's a full word match (ends line or followed by space)
                    if len(line) == name_len_in_line or (len(line) > name_len_in_line and line[name_len_in_line].isspace()):
                        # 3c. Extract sequence part
                        seq_part = line[name_len_in_line:].strip()
                        
                        # 3d. Append to taxon_to_sequence_parts
                        if actual_taxon_name_in_line not in taxon_to_sequence_parts:
                            taxon_to_sequence_parts[actual_taxon_name_in_line] = []
                        taxon_to_sequence_parts[actual_taxon_name_in_line].append(seq_part)
                        # found_taxon_for_line = True # Keep for debugging if needed
                        break # Found taxon for this line, move to next line
            # if not found_taxon_for_line: # Keep for debugging if needed
            #     print(f"Warning: Line orphaned or taxon name mismatch: {line[:50]}...")


        # 4. Assemble sequences from parts
        assembled_sequences = {}
        for taxon_original, parts in taxon_to_sequence_parts.items():
            full_sequence = "".join(parts)
            full_sequence = SEQUENCE_WHITESPACE_PATTERN.sub('', full_sequence) # Remove all whitespace
            assembled_sequences[taxon_original] = full_sequence

        # 5. Build final sequences list and 6. Handle missing taxa
        for taxon_original, taxon_clean in taxon_names:
            if taxon_original in assembled_sequences and assembled_sequences[taxon_original]:
                sequences.append((taxon_clean, assembled_sequences[taxon_original]))
            else:
                # This covers taxa in TAXLABELS but not in MATRIX or with empty sequence
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
        print(f"An error occurred during NEXUS parsing or FASTA writing: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print("Usage: python convert.treebase.nexus.to.fasta.py input_file.nexus output_file.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' does not exist or is not a file.")
        sys.exit(1)
    
    nexus_to_fasta(input_file, output_file)

if __name__ == "__main__":
    main()
