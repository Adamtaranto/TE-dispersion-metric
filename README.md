# TE-dispersion-metric

**Aim: **
To measure the degree to which repetitive elements are clustered or dispersed within a chromosome.
It has been widely observed that large genomic islands of repetative sequences (transposons etc) comprise 
isolate specific regions in fungal genomes and often contain virulence determinants. Similarly, many fungal
pathogens possess accessory chromosomes which are predominantly composed of repetative elements.  


Currently, transposon content is crudely summaraised as a proportion of the chromosome or genome space.
This measure is agnostic to the distribution of repetative sequences across the total genome space and is
therefore unable to differentiate between TE-poor genomes with a small number of repeat-islands or accessory 
chromosomes and genomes with a diffuse distribution of small repeats and lacking any isolated clusters.  


**Proposed metric:**

Within a chromosome - the mean distance from each repeat to all other repeats weighted by target feature's 
contribution to the total feature space.

```
## Chrom length
S = Total sequence length
## Distance between features scaled to chromosome length
d = distance between features / S
## Total space within S occupied by repeat features
T = Sum length of all transposons in S
## Weight of repeat feature as proportion of total repeat space T
w = feature length / T
## Count of non-overlapping features in Chrom
n = number of discrete features
```

Mean weighted distance from transposon i to all other features j

Mi = sum(d(i,j) * jw) / n - 1 

Mean weighted distance across all features i

sum(Mi) / n

## Example implementation

Requirements: 
  - biopython
  - pybedtools

```bash
./TEDcalc.py --gff repeats.gff3 --genome genome.fa
```
