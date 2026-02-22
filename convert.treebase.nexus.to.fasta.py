#!/usr/bin/env python3

# Convert a NEXUS file from TreeBASE into a FASTA file that lots of DNA analysis software can understand.

# Handles special characters well, at least the ones I have come across so far.   Older versions of this code didn't.

# By Alan Rockefeller - February 21, 2026

# Version 1.3


import sys
import re
import os

# Compiled regex patterns for use in nexus_to_fasta
# Support single and double quoted names, or unquoted tokens
TAXON_NAME_SPLIT_PATTERN = re.compile(r"'[^']*'|\"[^\"]*\"|\S+")
SEQUENCE_WHITESPACE_PATTERN = re.compile(r"\s+")


def unquote_taxon_name(name):
    """Removes single or double quotes from a taxon name if present."""
    if (name.startswith("'") and name.endswith("'")) or (
        name.startswith('"') and name.endswith('"')
    ):
        if len(name) > 1:
            return name[1:-1]
    return name


def strip_nexus_comments(text):
    """Removes bracketed NEXUS comments, handling nested brackets and quotes."""
    result = []
    depth = 0
    in_quote = None  # None, "'", or '"'

    for char in text:
        if in_quote:
            if char == in_quote:
                in_quote = None
            result.append(char)
        elif char == "'" or char == '"':
            in_quote = char
            result.append(char)
        elif char == '[':
            depth += 1
        elif char == ']':
            if depth > 0:
                depth -= 1
        elif depth == 0:
            result.append(char)
    return "".join(result)


def extract_nexus_block(content, block_name):
    """Extracts a block like MATRIX or TAXLABELS, handling nested brackets and semicolons inside quotes."""
    # Use word boundaries to avoid matching substrings of other words
    pattern = re.compile(r"\b" + re.escape(block_name) + r"\b\s*", re.IGNORECASE)
    match = pattern.search(content)
    if not match:
        return None
    
    start_index = match.end()
    depth = 0
    in_quote = None

    for i in range(start_index, len(content)):
        char = content[i]
        
        if in_quote:
            if char == in_quote:
                in_quote = None
        elif char == "'" or char == '"':
            in_quote = char
        elif char == '[':
            depth += 1
        elif char == ']':
            if depth > 0:
                depth -= 1
        elif char == ';' and depth == 0:
            return content[start_index:i].strip()
    return None


def make_unique(name, used_names):
    """Makes a FASTA header unique by appending _2, _3, etc. if needed."""
    if name not in used_names:
        used_names.add(name)
        return name
    base = name
    counter = 2
    while f"{base}_{counter}" in used_names:
        counter += 1
    new_name = f"{base}_{counter}"
    used_names.add(new_name)
    return new_name


def nexus_to_fasta(input_file, output_file):
    """
    Converts DNA sequences from NEXUS format to FASTA format.

    Args:
        input_file (str): Path to the input NEXUS file
        output_file (str): Path to the output FASTA file
    """
    content = None
    try:
        with open(input_file, "r", encoding="utf-8") as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

    if content is None:
        return

    # Find the MATRIX section in the NEXUS file
    matrix_content_raw = extract_nexus_block(content, "MATRIX")
    if matrix_content_raw is None:
        print("Error: Could not find the MATRIX section in the NEXUS file.")
        sys.exit(1)
        return

    # 1. Fix: Strip comments from the entire block at once to handle multi-line comments
    matrix_content_cleaned = strip_nexus_comments(matrix_content_raw)
    matrix_lines = matrix_content_cleaned.split("\n")

    try:
        # Find all taxon names in the TAXLABELS section
        taxlabels_content = extract_nexus_block(content, "TAXLABELS")
        taxon_names_info = []  # List of (original_name_for_matching, raw_name_from_file)
        has_taxlabels = False
        
        if taxlabels_content:
            has_taxlabels = True
            taxlabels_content = strip_nexus_comments(taxlabels_content)
            raw_names = TAXON_NAME_SPLIT_PATTERN.findall(taxlabels_content)
            for name_in_label in raw_names:
                original_name_for_matching = unquote_taxon_name(name_in_label)
                taxon_names_info.append((original_name_for_matching, name_in_label))

        taxon_to_sequence_parts = {}
        current_taxon_index_in_block = 0

        # Optimization: Sort taxa by length descending once to avoid re-sorting on every line
        sorted_taxa = sorted(taxon_names_info, key=lambda x: len(x[0]), reverse=True)
        # Optimization: Cache taxon names for index lookup
        taxon_list = [t[0] for t in taxon_names_info]

        for line in matrix_lines:
            line = line.strip()
            if not line:
                # 2. Fix: Blank lines reset the positional counter but don't stop discovery
                current_taxon_index_in_block = 0
                continue

            found_taxon = None
            name_len_in_line = 0

            token_match = TAXON_NAME_SPLIT_PATTERN.match(line)
            has_space_after_token = token_match and token_match.end() < len(line) and line[token_match.end()].isspace()

            # Match with existing known taxa using token-based approach
            if has_space_after_token:
                token = token_match.group(0)
                unquoted_token_cf = unquote_taxon_name(token).casefold()
                
                for taxon_original, _ in sorted_taxa:
                    if unquoted_token_cf == taxon_original.casefold():
                        found_taxon = taxon_original
                        name_len_in_line = token_match.end()
                        # Update positional tracker for interleaved consistency
                        current_taxon_index_in_block = taxon_list.index(found_taxon) + 1
                        break

            # Fallback: Discovery or Positional assignment
            if not found_taxon:
                # 1. Discovery (if allowed and looks like name+seq)
                if has_space_after_token and not has_taxlabels:
                    raw_name = token_match.group(0)
                    name_len_in_line = token_match.end()
                    found_taxon = unquote_taxon_name(raw_name)
                    if found_taxon not in taxon_list:
                        taxon_names_info.append((found_taxon, raw_name))
                        taxon_list.append(found_taxon)
                        # Re-sort only when new taxon is discovered
                        sorted_taxa = sorted(taxon_names_info, key=lambda x: len(x[0]), reverse=True)
                    
                    current_taxon_index_in_block = taxon_list.index(found_taxon) + 1
                
                # 2. Positional Assignment (Interleaved continuation)
                elif current_taxon_index_in_block < len(taxon_list):
                    # Check if the line is JUST a known taxon name (header-only line)
                    is_header_only = False
                    if token_match and token_match.end() == len(line):
                        token = token_match.group(0)
                        unquoted_token_cf = unquote_taxon_name(token).casefold()
                        for taxon_original in taxon_list:
                            if unquoted_token_cf == taxon_original.casefold():
                                # It's a header line: set index for next line, skip assignment
                                current_taxon_index_in_block = taxon_list.index(taxon_original)
                                is_header_only = True
                                break
                    
                    if not is_header_only:
                        # Interleaved positional assignment
                        found_taxon = taxon_list[current_taxon_index_in_block]
                        name_len_in_line = 0
                        current_taxon_index_in_block += 1
                
                # 3. Warning (looks like name+seq but no known taxon and no positional slot)
                elif has_space_after_token:
                    print(f"Warning: Line does not match any known taxon: {line[:50]}...")

            if found_taxon:
                seq_part = line[name_len_in_line:].strip()
                if found_taxon not in taxon_to_sequence_parts:
                    taxon_to_sequence_parts[found_taxon] = []
                taxon_to_sequence_parts[found_taxon].append(seq_part)

        # Assemble sequences and ensure unique FASTA headers
        used_fasta_names = set()
        final_sequences = []
        
        for taxon_original, raw_name in taxon_names_info:
            if taxon_original in taxon_to_sequence_parts:
                parts = taxon_to_sequence_parts[taxon_original]
                full_sequence = "".join(parts)
                full_sequence = SEQUENCE_WHITESPACE_PATTERN.sub("", full_sequence)
                
                clean_name = raw_name.replace("'", "").replace('"', "")
                clean_name = re.sub(r"\s+", "_", clean_name)
                clean_name = make_unique(clean_name, used_fasta_names)
                
                final_sequences.append((clean_name, full_sequence))
            else:
                print(f"Warning: No sequence found for taxon '{taxon_original}'")

        # Write sequences in FASTA format
        with open(output_file, "w", encoding="utf-8") as f:
            for taxon_name, sequence in final_sequences:
                f.write(f">{taxon_name}\n")
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i : i + 60] + "\n")

        print(f"Successfully converted {len(final_sequences)} sequences to FASTA format.")

    except Exception as e:
        print(f"An error occurred during NEXUS parsing or FASTA writing: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: python convert.treebase.nexus.to.fasta.py input_file.nexus output_file.fasta"
        )
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' does not exist or is not a file.")
        sys.exit(1)

    nexus_to_fasta(input_file, output_file)


if __name__ == "__main__":
    main()
