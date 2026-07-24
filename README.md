# strandshift

**Genomic sequence-model invariance auditor.**

Audit sequence classifiers across reverse complements and shifted windows before trusting variant-effect scores.

![strandshift cover](demo/cover.png)

![strandshift workbench](demo/dashboard.png)

## Why this exists

DNA foundation models are increasingly used for sequence classification and variant-effect prediction. A model can still change its answer when the same local biology is reverse-complemented or moved a few bases inside the input window. `strandshift` turns those transformations into a reproducible release gate.

## What ships

- Forward, reverse-complement, and signed window-shift transforms
- Transparent proxy model plus a stable JSON audit contract for external predictions
- CLI, JSON API, responsive local workbench, Docker image, tests, and CI
- No API keys and no sequence uploads

## Run it end to end

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -e .
strandshift demo
strandshift serve
```

Open <http://127.0.0.1:8090>. Analyze JSON with `strandshift analyze input.json`.

## Demo result

The committed sequence is evaluated under 18 strand/shift transformations. The deliberately brittle proxy model exposes a large prediction spread and fails a 0.05 invariance budget; the full per-transform trace makes the failure reproducible.

## Research basis

- [Nucleotide Transformer: building and evaluating robust foundation models for human genomics](https://www.nature.com/articles/s41592-024-02523-z)

## Scope

This audits transformation consistency. It does not establish biological validity, calibration, causal variant effects, or clinical utility.

## Test

```bash
python -m unittest discover -s tests -v
```

MIT licensed.
