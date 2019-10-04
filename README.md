# An exploration of cannabis genetic data

Cannabis research has been picking up a lot of momentum in the last few years. By now there's a ton (>1.5TB) of publicly available genetic data on cannabis at various places like the [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/sra/?term=cannabis+sativa) (NCBI). Ever wonder what's actually in it? I sure have! As far as I can tell, there isn't an exposé of it anywhere. This repository includes pre-generated variant calling format (VCF) files of sequencing data from Lynch _et al_, Sawler _et al_, and Phylos Biosciences, as well as some `R` code to explore and visualize the data.

# Requirements

You'll need `R` and the following dependencies:

- adegenet
- ape
- vcfR
- dartR
- poppr
- magrittr
- plotly
