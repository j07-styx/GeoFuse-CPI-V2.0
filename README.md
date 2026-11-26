# GeoFuse-CPI V2.0

**Dedicated research record for the frozen GeoFuse-CPI V2.0 candidate.** This repository contains methodology, architecture documentation, aggregate benchmark evidence, provenance and reproducibility boundaries. It is separate from the production website repository and legacy DynaCPI experimental repositories.

> **Status:** frozen research-release candidate; not independently validated for production replacement and not a deployable model package.

## Frozen artifact identity

- Checkpoint SHA-256: `b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4`
- Prediction archive SHA-256: `b8f935e2115ead2b83e247477118adedd38559acf368b80eabafc160240f0555`
- One frozen checkpoint is used for the headline V2.0 ledger. No V363/V499/v6/V48 scores are mixed into it.

## Architecture at a glance

- **Protein tower:** ESM-2 650M embeddings (1,280D) projected to 48D.
- **Compound tower:** 2,048-bit Morgan fingerprints (radius 2; chirality enabled) projected to 48D.
- **Feature fusion:** concatenation of protein, compound, elementwise product and absolute difference (192D).
- **Shared trunk:** 192 → 96 → 48, with task-specific scalar heads and frozen similarity blending.
- **Native structure expert:** separate 704D geometry residual on a legacy 423D v5 2D anchor; four source-trained folds. It is **not** end-to-end jointly trained with the V2.0 2D trunk.

The released 2D inference implementation also contains a 128D auxiliary input named `g`, zero-filled by the portable Predictor; its training-time provenance is not fully documented. This is explicitly flagged rather than guessed. **GenSPARC is a benchmark/task head, not a GenSPARC input token.**

See [Architecture](docs/architecture.md) and [Methodology](docs/methodology.md).

## Headline results

### Double-cold comparisons

| Dataset | V2.0 RMSE | Matched reproduced PerceiverCPI RMSE |
|---|---:|---:|
| DAVIS (n=1,339) | **0.6216** | 0.6828 |
| KIBA (n=3,442) | **0.7864** | 0.9267 |
| METZ (n=1,416) | **0.8797** | 1.0187 |

These are comparisons against **locally reproduced** PerceiverCPI outputs, not the original paper's headline table values.

### GenSPARC comparison

V2.0 has **20/24 numerical metric wins** against rounded published RMSE/Pearson point estimates across 12 cohorts. These are not paired significance tests and do not establish universal superiority. Four metric misses are listed in [the benchmark ledger](results/benchmark-summary.md).

### Native geometry

- BioLiP3 pooled retrospective n=70: native RMSE 1.3184 vs calibrated 2D RMSE 1.5760.
- OpenBind single-target n=482: native RMSE 1.7986 vs calibrated 2D RMSE 1.7917 (slightly worse).
- The prespecified OpenBind REAL-vs-SHUFFLED MAE interval crosses zero.

**Do not claim universal 3D superiority or independent multi-target external validation.** The external validation gate is blocked for lack of an eligible new cohort.

## Reproducibility

The original local V2.0 release verifier and its inference smoke check both passed on the original Windows/Python 3.11 research environment. This public repository contains aggregate records and a metadata checker, not model weights or raw datasets. A fresh clone alone cannot reproduce all model predictions.

Run the public-record check:

```bash
python scripts/verify_public_records.py
```

See [Reproducibility](docs/reproducibility.md), [Data provenance](docs/data-provenance.md) and [Limitations](docs/limitations-and-claim-boundaries.md).

## Repository boundaries

This repository contains V2.0 only. It does not contain the production v6 + V48 service, V363/V499 or later fusion experiments, website/frontend/backend deployment code, raw datasets, protected labels, checkpoint binaries or per-example predictions.

No open-source license is assigned pending ownership and third-party compatibility review. See [Release scope](docs/release-scope.md).
