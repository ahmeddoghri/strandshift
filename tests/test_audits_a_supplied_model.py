"""The audit has to fail brittle models *and* pass invariant ones."""
import unittest

from strandshift.core import DEMO, analyze, demo_model


def invariant(sequence: str, motif: str) -> float:
    """Best match anywhere in the window: no orientation or position bias."""
    if len(sequence) < len(motif):
        return 0.0
    return max(
        sum(a == b for a, b in zip(sequence[start:start + len(motif)], motif)) / len(motif)
        for start in range(len(sequence) - len(motif) + 1)
    )


def leftmost_only(sequence: str, motif: str) -> float:
    """Scores only the left edge, so any shift moves the answer."""
    return sum(a == b for a, b in zip(sequence[:len(motif)], motif)) / len(motif)


class SuppliedModelTest(unittest.TestCase):
    def test_an_invariant_model_passes_the_audit(self):
        # Specificity: an audit that fails everything is not an audit.
        report = analyze(DEMO, model=invariant)
        self.assertTrue(report["audit_passed"])
        self.assertEqual(report["prediction_range"], 0.0)

    def test_a_position_dependent_model_fails(self):
        report = analyze(DEMO, model=leftmost_only)
        self.assertFalse(report["audit_passed"])

    def test_the_bundled_demo_model_still_fails_with_its_published_numbers(self):
        report = analyze(DEMO)
        self.assertFalse(report["audit_passed"])
        self.assertEqual(report["prediction_range"], 0.3273)
        self.assertEqual(report["max_consensus_deviation"], 0.1818)

    def test_turning_off_the_position_bias_makes_the_demo_model_invariant(self):
        report = analyze({**DEMO, "position_bias": 0.0})
        self.assertTrue(report["audit_passed"])
        self.assertEqual(demo_model("ACGTCGATTA", "ACGTCGATTA", 0.0), 1.0)


if __name__ == "__main__":
    unittest.main()
