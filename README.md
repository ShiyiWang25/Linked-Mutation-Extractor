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
It also indicated that they are likely derived from a founder strain in the drug-naive sample containing the 25 specific SNVs.

![Image](https://github.com/ShiyiWang25/Toolkits_for_Nanopore_sequencing_data/blob/main/Figures/Application.png)
(A) HIV genomes in the virological-failure sample are aligned to the drug-naive consensus sequence and visualized using Integrative Genomics Viewer (IGV). Alignment reveals 25 positions across the HIV gag-pol region where mutations are selected and enriched in the virological-failure sample. Vertical lines in different colors highlight these 25 positions. Additionally, the DRM RT M184V is marked using a black arrowhead. RT: reverse transcriptase. 
(B) The nucleotide combinations at the 25 gag-pol positions in individual HIV genomes from the drug-naive sample (top) and the virological failure sample (bottom) are shown using waffle plots. In the waffle plots, each square represents one HIV genome. HIV genomes with the same nucleotide combination are colored in the same color, and HIV genomes with distinct nucleotide combinations are colored in different colors. The red arrow indicates a rare nucleotide combination from the drug-naive sample (top) enriched in HIV genomes in the virological failure sample (bottom). DRM: drug resistance mutation. PWH: people living with HIV. ART: antiretroviral therapy. 
