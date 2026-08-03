import unittest

from strandshift.core import DEMO, analyze, demo_model, reverse_complement


class ReverseComplementMotifTest(unittest.TestCase):
    """analyze() scored the reverse-complement strand's transforms against
    the ORIGINAL motif string, not its reverse complement. Reverse-
    complementing a sequence reverse-complements every substring inside it,
    including a motif occurrence -- so the same underlying biological
    signal appears as reverse_complement(motif) on the flipped strand, not
    as the original motif spelling. Scoring against the wrong reference
    checks whether the literal forward-strand motif happens to reappear by
    chance (it structurally can't), manufacturing a near-guaranteed
    "failure" on every reverse-complement transform rather than testing
    genuine strand invariance."""

    def test_rc_strand_recognizes_the_rc_motif_at_the_unshifted_baseline(self):
        sequence = DEMO["sequence"]
        motif = DEMO["motif"]
        rc_sequence = reverse_complement(sequence)
        rc_motif = reverse_complement(motif)
        # The RC'd sequence should score a perfect match against the RC'd
        # motif -- the same way the forward sequence scores a perfect match
        # against the forward motif -- because it's the identical biological
        # signal, just read from the other strand.
        forward_score = demo_model(sequence, motif, DEMO["position_bias"])
        rc_score = demo_model(rc_sequence, rc_motif, DEMO["position_bias"])
        self.assertAlmostEqual(forward_score, rc_score, places=6)
        self.assertEqual(round(rc_score, 4), 1.0)

    def test_rc_strand_predictions_now_use_the_rc_motif(self):
        result = analyze(DEMO)
        rc_predictions = [item for item in result["predictions"] if item["strand"] == "reverse_complement"]
        # With no window shift, the reverse-complement strand should
        # recognize its own (correctly oriented) motif just as cleanly as
        # the forward strand recognizes its baseline.
        unshifted = next(item for item in rc_predictions if item["shift"] == 0)
        self.assertEqual(unshifted["score"], 1.0)

    def test_fixed_instability_is_smaller_but_the_audit_still_fails(self):
        """The deliberately position-sensitive proxy model still genuinely
        fails the invariance budget once the motif-orientation bug is
        fixed -- the qualitative conclusion holds -- but the reported
        instability magnitude should be meaningfully smaller than the
        buggy version's, since half the transforms were being scored
        against the wrong reference."""
        result = analyze(DEMO)
        self.assertFalse(result["audit_passed"])
        self.assertLess(result["prediction_range"], 0.5)
        self.assertGreater(result["prediction_range"], 0.2)


if __name__ == "__main__":
    unittest.main()
