# ATOM3D-LBA prior structural research and the V2.0 boundary

**Status:** Historical research appendix. The results below are not GeoFuse-CPI V2.0 benchmark scores and are not included in the V2.0 headline results ledger.

## Why ATOM3D matters to this project

ATOM3D-LBA was an important part of GeoFuse's earlier structural research. It provided experimentally bound protein–ligand complexes for developing atom-pair descriptors, pocket controls, distance-shuffling controls, scaffold-cold evaluation and ligand-pose sensitivity tests. Those experiments materially informed the design and interpretation of the later native-structure work.

The official LBA 30%-sequence-identity split contains 3,507 training, 466 validation and 490 test complexes. LBA predicts pK = -log(K), where the underlying affinity is Ki or Kd depending on the source record. See the [ATOM3D LBA benchmark](https://www.atom3d.ai/lba.html), the [ATOM3D NeurIPS Datasets and Benchmarks paper](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/c45147dee729311ef5b5c3003946c48f-Abstract-round1.html), and the [Zenodo dataset record](https://doi.org/10.5281/zenodo.4914718).

## Relationship to the frozen V2.0 release

- The V2.0 2D headline ledger reports DAVIS, KIBA, METZ and GenSPARC.
- The V2.0 native expert is a separate component trained from the HiQBind TRAIN349 source and evaluated retrospectively on BioLiP3 and OpenBind.
- ATOM3D is not listed as a training or evaluation cohort in the frozen V2.0 results ledger. No evaluation of the frozen V2.0 checkpoint on ATOM3D has been verified.
- The ATOM3D experiments below belong to earlier exploratory research branches. Their scores were not merged into V2.0, and must not be score-stitched into its headline results.
- The portable V2.0 release does not contain a complete raw ATOM3D PDB/SDF-to-423D/704D feature-generation pipeline.

This is a version boundary, not a claim that the earlier ATOM3D work was unimportant. It is kept separate so readers can see what informed the research without mistaking an earlier model's results for V2.0 performance.

## Official validation and test exposure

The archived records show that the official 466-complex validation split was evaluated once by an earlier, separate structural model. The official 490-complex test split was subsequently scored once by a different, earlier 3,343D LightGBM experiment with three seeds. Neither result belongs to the V2.0 checkpoint.

### One-time official validation result

| Earlier model arm | RMSE |
|---|---:|
| Native physical geometry | 1.5810 |
| Shuffled geometry | 1.5881 |
| Pocket-only | 1.5760 |

The native-minus-shuffled 95% interval was [-0.0321, 0.0185], and native-minus-pocket was [-0.0213, 0.0300]. Both cross zero. The validation record explicitly says the official test had not yet been accessed at that point.

### One-time official test result

| Earlier 3,343D LightGBM arm | RMSE |
|---|---:|
| No-pocket 2D | 1.5783 |
| Pocket-only | 1.3661 |
| Shuffled geometry | 1.3702 |
| Native physical geometry | **1.3572** |

The native-minus-pocket 95% interval was [-0.03135, 0.01326]; native-minus-shuffled was [-0.03418, 0.00691]. Both cross zero. The native point estimate was lower, but the prespecified geometry-specific comparisons were inconclusive. The comparison against no-pocket 2D was stronger, but it does not isolate the contribution of correctly paired geometry from the added bound-pocket information.

**The official test is already consumed.** The final integrity audit records no_TEST_re_read_or_new_scoring = true, identifies the 490 scored records, and explicitly prohibits rereading or additional official-test scoring. We therefore do not reopen it or run the V2.0 checkpoint against it. Doing so would not create an untouched V2.0 evaluation; it would be a post hoc reuse of a test set whose results have already informed prior research.

## Recovered ATOM3D research milestones

All results in this section use the original ATOM3D TRAIN source or previously opened official split results. The training labels were used repeatedly during exploratory development. They are historical development evidence, not independent prospective validation.

| Branch | Evaluation design | Main result | Interpretation |
|---|---|---|---|
| V13 | 950-fit / 250-held-out original TRAIN; 62 protein groups | Scratch 3D RMSE 2.0129 vs 2D LightGBM 2.3222; equal-group macro RMSE 1.4147 vs 1.2846 | Pooled ranking reversed under group-macro evaluation; pose-pretrained model was worse (RMSE 2.3988). |
| V14 | 666-fit / 196-held-out original TRAIN; matched-capacity geometry ablation | Real 1.9919, shuffled 1.9924, contact-only 1.9927 RMSE; paired group intervals crossed zero | Fine-distance geometry did not show a meaningful advantage over shuffled/contact controls. |
| V16 | 707-fit / 163-held-out original TRAIN; 90 heuristic protein groups | 2D anchor 1.6727, real 1.6596, shuffled 1.6640, pocket-only 1.5987 RMSE | Pocket-only was best; all paired group-bootstrap intervals crossed zero. A separate frozen-weight test showed geometry-input sensitivity, not an affinity win. |
| V27 | 3,507 original TRAIN; four protein-component group folds; five matched LightGBM arms | Full native RMSE 1.5475, full shuffled 1.5683, pocket-only 1.5706 | Native beat pocket-only on group-macro RMSE, but native-versus-shuffled uncertainty crossed zero. |
| V28–V29 | Same trained V27 weights; geometry-only donor/shuffle interventions | Native pooled MSE was 3.79% lower than matched donor and 10.29% lower than shuffled | Evidence of learned pair-specific feature sensitivity; donor hybrids included out-of-range spatial values and were not guaranteed physical poses. |
| V30 | 2,961 RDKit-valid original TRAIN complexes; Murcko-scaffold-disjoint folds | Native MSE 1.8932 vs shuffled 1.9523 (3.03% lower); scaffold-group CI excluded zero | Promising retrospective scaffold-cold result, but protein groups still overlapped between fitting and validation. |
| V31 | 630-row stricter subset with both scaffold and heuristic protein-component absent from fitting | Native MSE 2.2558 vs shuffled 2.2839 (1.23% lower); paired group intervals crossed zero | The native-versus-shuffled advantage was not robust on the stricter dual-cold subset. Native also beat pocket-only by 3.95%, with exploratory intervals excluding zero. |
| V33 | ESM-2 35M plus Morgan; original TRAIN scaffold folds | ESM2+native MSE 1.7652 vs ESM2+shuffled 1.7890 (1.33% lower); group CI crossed zero | Pretrained sequence features improved the source baseline, but the native-specific geometry contrast remained inconclusive. On the stricter dual-cold subset, shuffled was numerically better than native. |
| V34–V36 | Rigid-body and near-native ligand-coordinate perturbations with steric screening | Pose sensitivity was measurable on source scaffold-cold folds; strict dual-cold intervals crossed zero | Useful engineering/falsification work, not evidence of reliable native-pose-specific affinity gains on unfamiliar proteins and scaffolds. |

### What these experiments established

1. The project developed several progressively stronger structural representations and meaningful pocket, shuffled-distance, donor and coordinate-perturbation controls.
2. Some native-geometry configurations improved over pocket-only or other baselines on particular retrospective source splits.
3. The strongest geometry-specific conclusions weakened under stricter scaffold-plus-protein holdout, and several paired confidence intervals crossed zero.
4. The V33 sequence experiment showed that stronger pretrained protein features improved the source development baseline, but did not establish a robust native-versus-shuffled gain.
5. V34–V36 demonstrated model sensitivity to ligand-coordinate changes. Sensitivity alone does not establish that the correct pose improves affinity prediction on new targets.

The full experimental sequence, protocols, result files and local evidence hashes are recorded in the source manifest below. Raw ATOM3D structures, labels, per-example predictions and derived feature arrays are not redistributed here.

## Scientific conclusion

ATOM3D was an important structural-development dataset for GeoFuse. It helped build and falsify geometry representations, quantify the role of pocket information, test scaffold generalization and expose the difference between geometric input sensitivity and a genuine affinity-prediction benefit.

The defensible conclusion is **mixed retrospective evidence of structural utility**, not universal 3D superiority. The one-time official test point estimate favored native geometry, but its native-versus-pocket and native-versus-shuffled intervals crossed zero. It was generated by an earlier LightGBM model, not V2.0.

No ATOM3D score is added to the V2.0 headline benchmark ledger. A new independent evaluation requires a genuinely new, unopened, assay-verified multi-target cohort with source-lineage, scaffold and protein-family audits, physically credible pose controls, and frozen models before labels are accessed.

## Source provenance and licensing

The official ATOM3D LBA dataset is distributed under a Creative Commons Attribution-NonCommercial-NoDerivatives (CC BY-NC-ND) license. This repository does not include the original ATOM3D archives, structures, numerical labels, label-bearing arrays, derived feature arrays or per-example predictions. Obtain data directly from the [official Zenodo record](https://doi.org/10.5281/zenodo.4914718) and comply with its current license and citation requirements.

The local research evidence hashes are recorded in [manifests/atom3d-evidence-manifest.json](../manifests/atom3d-evidence-manifest.json). The compact aggregate record is [results/atom3d-historical-summary.json](../results/atom3d-historical-summary.json). These records establish provenance and internal consistency; they do not rerun training or reconstruct the original benchmark from raw data.

Run python scripts/verify_atom3d_summary.py from the repository root to check the aggregate summary and provenance flags. The verifier reads only the two compact JSON records; it does not access ATOM3D data or labels and does not rescore the official test.
