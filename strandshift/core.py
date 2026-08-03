from __future__ import annotations

from statistics import mean, pstdev


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


DEMO = {
    "sequence": "NNNNNNNNNNNACGTCGATTANNNNNNNNNNN",
    "motif": "ACGTCGATTA",
    "max_shift": 4,
    "position_bias": 0.9,
}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def shift_window(sequence: str, amount: int) -> str:
    if amount > 0:
        return ("N" * amount + sequence)[: len(sequence)]
    if amount < 0:
        return (sequence[-amount:] + "N" * (-amount))[: len(sequence)]
    return sequence


def _similarity(left: str, right: str) -> float:
    return sum(a == b for a, b in zip(left, right)) / max(1, len(right))


def demo_model(sequence: str, motif: str, position_bias: float) -> float:
    """A transparent proxy for a strand- and position-sensitive sequence model."""
    if len(sequence) < len(motif):
        return 0.0
    center = (len(sequence) - len(motif)) / 2
    scores = []
    for start in range(len(sequence) - len(motif) + 1):
        match = _similarity(sequence[start : start + len(motif)], motif)
        distance = abs(start - center) / max(1.0, center)
        scores.append(match * (1.0 - position_bias * min(1.0, distance)))
    return max(scores)


def analyze(payload: dict) -> dict:
    sequence = str(payload["sequence"]).upper()
    motif = str(payload["motif"]).upper()
    max_shift = int(payload.get("max_shift", 3))
    position_bias = float(payload.get("position_bias", 0.5))
    if not sequence or set(sequence) - set("ACGTN"):
        raise ValueError("sequence must contain only A, C, G, T, or N")
    if not motif or set(motif) - set("ACGT"):
        raise ValueError("motif must contain only A, C, G, or T")
    if not 0 <= max_shift <= 20:
        raise ValueError("max_shift must be between 0 and 20")
    if not 0 <= position_bias <= 1:
        raise ValueError("position_bias must be between 0 and 1")

    predictions = []
    # On the reverse-complement strand, the SAME underlying biological motif
    # appears as its reverse complement, not as the original motif string --
    # reverse-complementing a sequence reverse-complements every substring
    # inside it, including the motif occurrence. Scoring the RC'd sequence
    # against the un-RC'd motif checks whether the literal forward-strand
    # spelling happens to reappear (it structurally can't, except by chance),
    # which isn't a test of strand invariance at all -- it manufactures a
    # near-guaranteed "failure" on every reverse-complement transform.
    for strand, oriented, oriented_motif in (
        ("forward", sequence, motif),
        ("reverse_complement", reverse_complement(sequence), reverse_complement(motif)),
    ):
        for shift in range(-max_shift, max_shift + 1):
            transformed = shift_window(oriented, shift)
            predictions.append(
                {
                    "strand": strand,
                    "shift": shift,
                    "score": round(demo_model(transformed, oriented_motif, position_bias), 4),
                }
            )

    scores = [item["score"] for item in predictions]
    baseline = next(item["score"] for item in predictions if item["strand"] == "forward" and item["shift"] == 0)
    consensus = mean(scores)
    max_deviation = max(abs(score - consensus) for score in scores)
    return {
        "sequence_length": len(sequence),
        "transform_count": len(predictions),
        "baseline_score": round(baseline, 4),
        "invariant_consensus": round(consensus, 4),
        "prediction_range": round(max(scores) - min(scores), 4),
        "prediction_std": round(pstdev(scores), 4),
        "max_consensus_deviation": round(max_deviation, 4),
        "audit_passed": max_deviation <= 0.05,
        "worst_transform": max(predictions, key=lambda item: abs(item["score"] - baseline)),
        "predictions": predictions,
        "scope": "Transformation audit for sequence classifiers; biological validity and calibration require task-specific data.",
    }
