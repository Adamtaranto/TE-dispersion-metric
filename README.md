# TE-dispersion-metric
Proposed metric for measuring the clustering of transposons within a genome

Inverse mean distance to all other features weighted by target feature contribution to total feature space.

S = Total sequence length
d = distance between features / S
T = Sum length of all transposons in S
w = feature length / T
n = number of discrete features

Mean distance from each transposon i to all other features j

Mi = d(i,j) * jw / n - 1 

Mean weighted distance across all features

∑ Mi / n

## Example implementation

Requirements: 
  - biopython
  - pybedtools

```bash
./TEDcalc.py --gff repeats.gff3 --genome genome.fa
```
