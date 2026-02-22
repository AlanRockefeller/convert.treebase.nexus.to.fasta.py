import unittest
import os
import sys
import io
import tempfile
import shutil
from unittest.mock import patch
import importlib.util

# Dynamically import the module
# Assuming test_converter.py is in the root of the repository,
# and convert.treebase.nexus.to.fasta.py is also in the root.
module_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "convert.treebase.nexus.to.fasta.py")
)
module_name = "converter_module"

spec = importlib.util.spec_from_file_location(module_name, module_path)
converter_module = importlib.util.module_from_spec(spec)
# It's important to add the module to sys.modules BEFORE exec_module,
# especially if the module itself has internal imports that might rely on its own name.
sys.modules[module_name] = converter_module
spec.loader.exec_module(converter_module)

nexus_to_fasta = converter_module.nexus_to_fasta
converter_main = converter_module.main


class TestNexusToFastaConverter(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.input_file_path = os.path.join(self.test_dir, "test_input.nexus")
        self.output_file_path = os.path.join(self.test_dir, "test_output.fasta")

        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        sys.stdout = io.StringIO()
        sys.stderr = io.StringIO()
        # No need to flush StringIO initially

    def tearDown(self):
        sys.stdout.flush()  # Flush before restoring
        sys.stderr.flush()
        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr
        shutil.rmtree(self.test_dir)

    def _write_input_nexus(self, content):
        with open(self.input_file_path, "w", encoding="utf-8") as f:
            f.write(content)

    def _read_output_fasta(self):
        if os.path.exists(self.output_file_path):
            with open(self.output_file_path, "r", encoding="utf-8") as f:
                return f.read()
        return None

    def test_simple_conversion(self):
        nexus_content = """
#NEXUS
BEGIN TAXA;
    TAXLABELS T1 T2 T3;
END;
BEGIN CHARACTERS;
    MATRIX
    T1 ACGT
    T2 TGCA
    T3 AAAA
    ;
END;
"""
        self._write_input_nexus(nexus_content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)

        self.assertTrue(os.path.exists(self.output_file_path))
        fasta_content = self._read_output_fasta()
        expected_fasta = ">T1\nACGT\n>T2\nTGCA\n>T3\nAAAA\n"
        self.assertEqual(fasta_content, expected_fasta)

    def test_long_sequence_wrapping(self):
        seq_part1 = "A" * 60
        seq_part2 = "C" * 20
        nexus_content = f"""
#NEXUS
BEGIN TAXA;
    TAXLABELS LongSeqTaxon;
END;
BEGIN CHARACTERS;
    MATRIX
    LongSeqTaxon {seq_part1}{seq_part2}
    ;
END;
"""
        self._write_input_nexus(nexus_content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        expected_fasta = f">LongSeqTaxon\n{seq_part1}\n{seq_part2}\n"
        self.assertEqual(fasta_content, expected_fasta)

    def test_quoted_and_spaced_taxon_names(self):
        # Simplified to focus on the problematic "taxon two"
        nexus_content = """
#NEXUS
BEGIN TAXA; TAXLABELS "taxon two"; END;
BEGIN CHARACTERS; MATRIX "taxon two" TGCA; END;
"""
        self._write_input_nexus(nexus_content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)

        fasta_content = self._read_output_fasta()
        expected_fasta = ">taxon_two\nTGCA\n"  # Expected output for the simplified case
        self.assertEqual(fasta_content, expected_fasta)

    @patch(f"{module_name}.sys.exit", side_effect=SystemExit(1))  # Patch sys.exit inside the loaded module
    def test_missing_matrix_section(self, mock_sys_exit):
        nexus_content = """#NEXUS SIMPLIFIED CONTENT NOMATRIX"""  # Minimal content, no MATRIX, no semicolon
        self._write_input_nexus(nexus_content)

        with self.assertRaises(SystemExit):
            nexus_to_fasta(self.input_file_path, self.output_file_path)

        mock_sys_exit.assert_called_once_with(1)
        self.assertIn(
            "Error: Could not find the MATRIX section in the NEXUS file.",
            sys.stdout.getvalue(),
        )

    @patch(f"{module_name}.sys.exit", side_effect=SystemExit(1))  # Patch sys.exit inside the loaded module
    def test_input_file_not_found_main(self, mock_sys_exit):
        non_existent_file = os.path.join(self.test_dir, "no_such_file.nexus")

        original_argv = sys.argv
        sys.argv = ["converter_script_name", non_existent_file, self.output_file_path]
        try:
            with self.assertRaises(SystemExit):
                converter_main()
        finally:
            sys.argv = original_argv  # Ensure sys.argv is restored

        mock_sys_exit.assert_any_call(
            1
        )  # main() might call exit, then nexus_to_fasta might too
        expected_error_msg = (
            f"Error: Input file '{non_existent_file}' does not exist or is not a file."
        )
        # Check if the specific error from main() is present.
        all_stdout = sys.stdout.getvalue()  # Changed to stdout
        self.assertIn(
            expected_error_msg, all_stdout
        )  # Check specific message from main

    @patch(f"{module_name}.sys.exit", side_effect=SystemExit(1))
    # Removed @patch for print, will rely on sys.stdout redirection.
    def test_input_file_not_found_direct_call(
        self, mock_sys_exit
    ):  # mock_print removed
        non_existent_file = os.path.join(self.test_dir, "no_such_file.nexus")
        with self.assertRaises(SystemExit):
            nexus_to_fasta(non_existent_file, self.output_file_path)
        mock_sys_exit.assert_called_once_with(1)
        self.assertIn(
            f"Error: Input file not found: {non_existent_file}", sys.stdout.getvalue()
        )  # Changed to stdout

    def test_taxon_in_labels_not_in_matrix(self):
        nexus_content = """
#NEXUS
BEGIN TAXA;
    TAXLABELS T1 T2_missing T3;
END;
BEGIN CHARACTERS;
    MATRIX
    T1 ACGT
    T3 AAAA
    ;
END;
"""
        self._write_input_nexus(nexus_content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)

        fasta_content = self._read_output_fasta()
        expected_fasta = ">T1\nACGT\n>T3\nAAAA\n"  # T2_missing should not be in output
        self.assertEqual(fasta_content, expected_fasta)
        # Check stdout for the warning
        self.assertIn(
            "Warning: No sequence found for taxon 'T2_missing'", sys.stdout.getvalue()
        )

    def test_empty_matrix_no_taxa_in_matrix(self):
        nexus_content = """
#NEXUS
BEGIN TAXA;
    TAXLABELS T1 T2;
END;
BEGIN CHARACTERS;
    MATRIX
    ;
END;
"""
        self._write_input_nexus(nexus_content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        self.assertEqual(
            fasta_content,
            "",
            "FASTA file should be empty if no sequences are processed.",
        )
        # Check stdout for warnings for both taxa
        self.assertIn(
            "Warning: No sequence found for taxon 'T1'", sys.stdout.getvalue()
        )
        self.assertIn(
            "Warning: No sequence found for taxon 'T2'", sys.stdout.getvalue()
        )

    def test_curly_brace_crash(self):
        # Issue 1: Curly brace in taxon name causes format() to crash
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS {A} B;
END;
BEGIN CHARACTERS;
    MATRIX
    {A} ACGT
    B   TGCA
    ;
END;
"""
        self._write_input_nexus(content)
        # Should not crash
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        self.assertIn(">{A}\nACGT\n", fasta_content)

    def test_quoted_comments(self):
        # Issue 2: Brackets inside quotes should not be treated as comments
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS "name[with]bracket";
END;
BEGIN CHARACTERS;
    MATRIX
    "name[with]bracket" ACGT
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        # The parser cleans quotes, so "name[with]bracket" -> name[with]bracket
        self.assertIn(">name[with]bracket\nACGT", fasta_content)

    def test_interleaved_continuation_with_taxlabels(self):
        # Issue 3: Interleaved lines without names should work even if TAXLABELS exists
        # Updated to include a blank line to reset the positional counter (new logic)
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS A B;
END;
BEGIN CHARACTERS;
    MATRIX
    A ACG
    B ACG

    T
    T
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        self.assertIn(">A\nACGT\n", fasta_content)
        self.assertIn(">B\nACGT\n", fasta_content)

    def test_case_insensitive_matching(self):
        # Issue 5: Case mismatch should still match
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS TaxonA;
END;
BEGIN CHARACTERS;
    MATRIX
    taxona ACGT
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        self.assertIn(">TaxonA\nACGT\n", fasta_content)

    def test_no_wrap_around_within_block(self):
        # Verify that positional fallback does NOT wrap around within a single block
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS A B;
END;
BEGIN CHARACTERS;
    MATRIX
    A ACG
    B ACG
    T
    T
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        # Since there's no blank line, the two 'T' lines should be ignored
        self.assertIn(">A\nACG\n", fasta_content)
        self.assertIn(">B\nACG\n", fasta_content)
        self.assertNotIn("ACGT", fasta_content)

    def test_casefold_matching(self):
        # Issue: .casefold() should handle Unicode characters more robustly than .lower()
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS Taxonß;
END;
BEGIN CHARACTERS;
    MATRIX
    TaxonSS ACGT
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        # Taxonß.casefold() == taxonss
        self.assertIn(">Taxonß\nACGT\n", fasta_content)

    def test_name_on_separate_line(self):
        # Support for taxon name on its own line followed by sequence
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS A B;
END;
BEGIN CHARACTERS;
    MATRIX
    A
    ACG
    B
    TGC
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        self.assertIn(">A\nACG\n", fasta_content)
        self.assertIn(">B\nTGC\n", fasta_content)

    def test_no_empty_chunks_from_names(self):
        # Ensure that name-only lines don't add empty strings to sequence parts
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS A;
END;
BEGIN CHARACTERS;
    MATRIX
    A
    ACG
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        # If an empty chunk was added, it might affect sequence length or formatting
        self.assertEqual(fasta_content, ">A\nACG\n")

    def test_interleaved_spaced_sequence(self):
        # Support for interleaved continuation lines that contain spaces
        content = """#NEXUS
BEGIN TAXA;
    TAXLABELS A B;
END;
BEGIN CHARACTERS;
    MATRIX
    A ACG
    B ACG

    A T G C
    B T G C
    ;
END;
"""
        self._write_input_nexus(content)
        nexus_to_fasta(self.input_file_path, self.output_file_path)
        fasta_content = self._read_output_fasta()
        # Spaces should be stripped during assembly
        self.assertIn(">A\nACGTGC\n", fasta_content)
        self.assertIn(">B\nACGTGC\n", fasta_content)


if __name__ == "__main__":
    unittest.main()
