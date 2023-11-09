# Toolkits_for_Nanopore_sequencing_data
This repo contains Python-based analysis tools for analyzing long-read sequencing data.
---

## Linked-Mutation Extractor (LiME)
## Description
This _Python_ algorithm performs concurrent variant calling across multiple loci in aligned sequences and clusters the sequences based on their unique nucleotide combinations (genetic patterns).
The output is visualized using a Waffle plot, wherein each square represents a sequence, and a cluster of sequences characterized by the same nucleotide combination is colored the same.

## Features
1. LiME takes the BAM file as input.
2. LiME generates a waffle plot using the [PyWaffle](https://pywaffle.readthedocs.io/en/latest/) package.
3. LiME stores the analyzing outcome in a CSV file, which contains detailed information about genetic patterns and the number and ID of the reads harboring each genetic pattern.

## Application
This algorithm has been utilized in [(Gallardo et al. 2021)](https://academic.oup.com/nar/article/49/12/e70/6225234) to illustrate potential mechanisms that gave rise to the enrichment for 25 specific SNVs within the HIV genome during virological failure.
It shows that the 25 SNVs are likely collectively enriched in the same HIV genome instead of being selected separately in different genomes.
