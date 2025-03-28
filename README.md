# TreeBASE NEXUS to FASTA Converter

A Python utility for converting DNA sequence alignments downloaded from TreeBASE in NEXUS format to FASTA format.

*By Alan Rockefeller - March 28, 2025*

## Overview

This tool converts NEXUS files downloaded from TreeBASE (and other NEXUS files too) into FASTA format, which is widely supported by DNA analysis software. It handles special characters and quoted taxon names that can often cause problems in conversion.

## Features

- Properly extracts DNA sequences from NEXUS formatted files downloaded from TreeBASE
- Correctly handles taxon names with spaces or special characters
- Preserves quoted names as single entities
- Formats FASTA sequences in standard 60-character blocks
- Produces clean FASTA headers that are compatible with most bioinformatics software

## Usage

```bash
python convert.treebase.nexus.to.fasta.py input_file.nexus output_file.fasta
```

## How It Works

1. The script first extracts all taxon names from the TAXLABELS section of the NEXUS file
2. It creates clean versions of these names for FASTA (replacing quotes and spaces with underscores)
3. It then finds each sequence in the matrix by looking for the original taxon name
4. Finally, it writes each sequence in standard FASTA format with 60 characters per line


## Common Issues and Solutions

If you encounter problems with converted FASTA files in MEGA or other software:

1. **Access violation errors**: These can occur when building a phylogenetic tree in MEGA 12 when sequence headers contain special characters, causing the names to spill over into the nucleotides. The script handles this by cleaning taxon names.

2. **Ambiguity codes**: While the script preserves standard IUPAC ambiguity codes (R, Y, M, K, S, W, etc.), some software may have issues with these. Ambiguity codes contain important information, so if your bioinformatics software doesn't handle them, it would be best to find software that does.

3. **Multi-line sequences**: The script formats sequences in standard 60-character blocks for optimal compatibility.

## Requirements

- Python 3.x
- No additional libraries required

## License

This code is provided under the MIT License.
