import unittest

from strandshift.core import DEMO, analyze, reverse_complement, shift_window


class StrandShiftTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AACGTN"), "NACGTT")

    def test_shift_preserves_window_length(self):
        self.assertEqual(shift_window("ACGT", 2), "NNAC")
        self.assertEqual(shift_window("ACGT", -2), "GTNN")

    def test_demo_surfaces_instability(self):
        result = analyze(DEMO)
        self.assertEqual(result["transform_count"], 18)
        self.assertGreater(result["prediction_range"], 0.25)
        self.assertFalse(result["audit_passed"])

    def test_invalid_sequence_fails(self):
        with self.assertRaisesRegex(ValueError, "sequence"):
            analyze({**DEMO, "sequence": "ACGTX"})


if __name__ == "__main__":
    unittest.main()
