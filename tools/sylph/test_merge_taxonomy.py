"""Regression tests for zero-hit samples and real upstream merge failures.

Run with sylph-tax 1.2.0 on PATH: python -m unittest test_merge_taxonomy.py
"""

import subprocess
import tempfile
import unittest
from pathlib import Path

from merge_taxonomy import merge_taxonomy


class MergeTaxonomyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.profile = self.root / "profile.tsv"
        self.profile.write_text("Sample_file\tGenome_file\n")
        self.sketches = self.root / "sketches"
        self.sketches.mkdir()
        self.taxonomy = self.root / "taxonomy"
        self.taxonomy.mkdir()
        self.output = self.root / "merged.tsv"

    def merge(self, paired=False):
        merge_taxonomy(self.profile, self.sketches, self.taxonomy, self.output, "ANI", paired)

    def positive_sample(self):
        self.profile.write_text("Sample_file\tGenome_file\npositive.fastq.gz\tgenome.fasta\n")
        (self.sketches / "positive.fastq.gz.sylsp").touch()
        (self.taxonomy / "positive.sylphmpa").write_text(
            "#SampleID\tpositive.fastq.gz\n"
            "clade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)\n"
            "d__Bacteria|t__genome\t100.0\t100.0\t99.95\t0.047\n"
        )

    def test_all_negative_samples_retain_columns(self):
        for name in ["SRR24951297.fastq.gz", "SRR24951347.fastq.gz"]:
            (self.sketches / f"{name}.sylsp").touch()
        self.merge()
        self.assertEqual(
            self.output.read_text(),
            "clade_name\tSRR24951297.fastq.gz\tSRR24951347.fastq.gz\n",
        )

    def test_negative_paired_sample_name(self):
        (self.sketches / "SRR24951297.fastq.gz.paired.sylsp").touch()
        self.merge(paired=True)
        self.assertEqual(self.output.read_text(), "clade_name\tSRR24951297.fastq.gz\n")

    def test_mixed_samples_retain_positive_values_and_negative_column(self):
        self.positive_sample()
        (self.sketches / "negative.fastq.gz.sylsp").touch()
        self.merge()
        self.assertEqual(
            self.output.read_text(),
            "clade_name\tnegative.fastq.gz\tpositive.fastq.gz\n"
            "d__Bacteria|t__genome\t0\t99.95\n",
        )

    def test_empty_profile_is_not_a_biological_negative(self):
        (self.sketches / "sample.fastq.gz.sylsp").touch()
        self.profile.write_text("")
        with self.assertRaisesRegex(ValueError, "Missing Sylph profile header"):
            self.merge()
        self.assertFalse(self.output.exists())

    def test_upstream_merge_failure_propagates(self):
        self.positive_sample()
        (self.taxonomy / "positive.sylphmpa").write_text("malformed taxonomy\n")
        with self.assertRaises(subprocess.CalledProcessError):
            self.merge()


if __name__ == "__main__":
    unittest.main()
