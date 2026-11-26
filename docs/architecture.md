# GeoFuse-CPI V2.0 architecture

## Scope and source of truth

This page describes the frozen V2.0 research-release implementation. It does not describe V363/V499 experiments, the production v6 service, or a proposed V2.0+V363 fusion.

- Checkpoint SHA-256: `b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4`
- Frozen prediction archive SHA-256: `b8f935e2115ead2b83e247477118adedd38559acf368b80eabafc160240f0555`

## 1. 2D sequence-and-chemistry predictor

### Inputs

- Protein: ESM-2 650M embedding, 1,280 dimensions.
- Compound: Morgan fingerprint, 2,048 bits, radius 2, chirality enabled.

### Protein tower

The released implementation applies a 1,280 → 128 → 48 projection with LayerNorm, GELU and dropout in the hidden projection.

There is also an auxiliary 128-dimensional input named `g`. Its projection is 128 → 48 → 48 and is added to the protein representation. The released portable Predictor currently supplies an all-zero `g` vector at inference. The public release notes do not fully specify the training-time provenance/semantics of this auxiliary input. It must not be relabelled as a GenSPARC token or as native geometry without additional source evidence.

### Compound tower

The Morgan fingerprint is projected 2,048 → 128 → 48, with LayerNorm, GELU and dropout in the hidden projection.

### Interaction feature fusion

Let (p) be the 48-dimensional protein-side representation and (c) the 48-dimensional compound representation. The released implementation concatenates:

- (p) (48)
- (c) (48)
- elementwise product (p \odot c) (48)
- absolute difference (|p-c|) (48)

This produces a **192-dimensional** interaction vector.

### Shared trunk and task heads

The shared MLP is 192 → 96 → 48, followed by assay/task-specific scalar heads for DAVIS, KIBA, METZ and GenSPARC. The public raw-input Predictor explicitly exposes DAVIS, KIBA and METZ inference; the GenSPARC head is part of the benchmark-task implementation, not a special input token.

The final assay prediction also blends the neural output with a task-specific similarity estimate using a frozen kernel weight. The similarity estimate uses Morgan Tanimoto similarity and centered protein-embedding similarity with frozen training anchors. DAVIS additionally applies the documented output floor.

## 2. Native-structure expert

The native expert is a **separate structural residual component**, not a jointly trained continuation of the V2.0 48-dimensional protein/compound trunk.

Its inputs are:
- a 423-dimensional legacy 2D context: 422 dimensions of frozen v5 tower context plus the raw v5 affinity prediction;
- a 704-dimensional native ligand–pocket geometry descriptor.

The geometry descriptor represents ligand–pocket atom-pair chemistry, radial-distance bases and bond/contact angular bases from experimentally observed bound structures.

For each of four source-trained folds, the native expert:
1. obtains a calibrated 2D base prediction from the fold's LightGBM calibrators;
2. standardizes the 704-dimensional geometry features;
3. applies a residual MLP 704 → 64 → 32 → 1, with LayerNorm and SiLU;
4. bounds the residual with `0.5 * tanh(raw_residual)`;
5. adds the scaled residual to the calibrated 2D base.

The four fold outputs are averaged. Native, shuffled-geometry and pocket-only arms are retained as distinct controls.

## 3. Accurate architecture description

Use: **a frozen hybrid research release containing a 2D protein/compound interaction predictor and a separately trained native-structure residual expert.**

Do not describe it as an end-to-end jointly trained 2D+3D network. Do not depict GenSPARC as an input token. GenSPARC is a benchmark/task in the released evaluation; it is not the native ligand–pocket geometry input.

## 4. Known implementation-documentation gap

The auxiliary `g` input is present in the released inference class, projected into the protein-side representation, and zero-filled by the portable Predictor. Its training-time feature provenance is not sufficiently described in the public release files. This is an explicit documentation/reproducibility item to resolve before making claims about that feature. It is intentionally not guessed here.
