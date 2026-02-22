## Changelog

### Version 1.3 (2-21-2026)

Improved NEXUS parsing robustness and compatibility, especially for TreeBASE exports and edge cases.

### Added
- **Nested comment handling** for NEXUS bracket comments (`[ ... ]`) via a comment-stripping parser.
- **Robust block extraction** for `MATRIX` and `TAXLABELS` that ignores semicolons inside comments.
- **Fallback parsing when `TAXLABELS` is missing** by discovering taxa from the first `MATRIX` block.
- **Duplicate FASTA header protection** with deterministic suffixes (`_2`, `_3`, etc.) when cleaned names collide.

### Changed
- `MATRIX` / `TAXLABELS` detection is now **case-insensitive** (works with `Matrix`, `matrix`, etc.).
- Taxon token parsing now supports:
  - single-quoted names (`'Taxon name'`)
  - double-quoted names (`"Taxon name"`)
  - unquoted names
- FASTA header normalization now replaces **all whitespace** with underscores (not just spaces).
- Interleaved matrix parsing behavior improved:
  - fallback taxon discovery only happens in the **first matrix block**
  - later unmatched continuation lines are no longer misread as taxa

### Fixed
- Parser no longer fails silently (0 sequences) when valid `MATRIX` data exists but `TAXLABELS` is absent.
- Double-quoted taxon names in `TAXLABELS` are now parsed correctly.
- Semicolons inside comments no longer prematurely terminate `MATRIX` / `TAXLABELS` extraction.
- Inline NEXUS comments no longer leak into FASTA sequence output.
- Duplicate cleaned taxon names no longer produce duplicate FASTA headers.

### Notes
- Output FASTA formatting is unchanged (60-character line wrapping preserved).
- CLI usage and general error/success messaging remain the same.
