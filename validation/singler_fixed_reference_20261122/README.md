# SingleR comparison with a fixed reference

This independent rerun changes only the marker-selection setting. No manuscript or figure was changed.

## Configuration

- Query: the same 4,913 GSE268919 cells in the saved sampling manifest, with 26 frozen cluster assignments.
- Evaluation: the same 10 locally reviewed clusters used for Figure 2. The local review is an exploratory reference, not an independent subtype ground truth.
- Both arms: SingleR 2.6.0, celldex 1.14.0 MouseRNAseqData, 358 reference profiles and 28 fine labels; cluster mode and pruned predictions.
- Default arm: automatic marker selection across shared genes.
- Dictionary arm: automatic marker selection restricted to 713 eligible PCOS-CAM marker genes via `restrict`, after one-to-one mouse ortholog mapping and intersection with query/reference genes.
- The delivered File_1.xlsx was used. Its parsed marker and ortholog tables were verified as identical to the frozen algorithm dictionary. Dictionary vote weights are not used by SingleR.
- Both arms retain the same reference expression profiles and reference-derived output labels. This comparison does not add ovarian labels to MouseRNAseqData.

## Results

| Arm | Reviewed clusters | Broad-family matches | Distinct finer labels among family-matched clusters |
|---|---:|---:|---:|
| SingleR with PCOS-CAM marker restriction | 10 | 6 | 4 |
| SingleR with automatic markers | 10 | 6 | 4 |

The same six clusters match the reviewed broad families: 2, 9, 11, 13, 16 and 20. Six reviewed clusters change their labels, but none changes broad-family match status. Stromal clusters 2 and 13 change from Fibroblasts to Fibroblasts activated. The four counted labels are B cells, Fibroblasts/Fibroblasts activated, Macrophages and NK cells. Endothelial cells is a family-only call and is excluded from the distinct finer-label count under the existing evaluation rules.

There are 14 changed labels across all 26 clusters. The default run exactly reproduces all 26 saved native predictions. Marker restriction does not improve either reported summary metric in this configuration. The finer-label count describes output resolution and does not validate subtype accuracy.

## Files

- `metrics.csv`: primary count summary.
- `reviewed_cluster_comparison.csv`: the 10 reviewed clusters with both predictions and counting flags.
- `cluster_predictions.csv`: predictions for all 26 clusters.
- `changed_clusters.csv`: the 14 clusters with changed labels.
- `counted_fine_labels.csv`: counted labels and cluster frequencies.
- `marker_genes.csv`: the 713 eligible mouse genes and human orthologs.
- `reference_sample_labels.csv`: labels for the same 358 reference profiles.
- `SingleR_default_scores.rds`, `SingleR_dictionary_scores.rds`: full SingleR outputs.
- `input_manifest.csv`, `run_configuration.txt`, `verification.txt`, `sessionInfo.txt`, `run.log`: inputs, checks and environment.

Run from the PCOS-CAM root:

```sh
Rscript validation/singler_fixed_reference_20261122/compare.R validation/subtype_evidence/query_input.rds
```
