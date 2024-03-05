# RegionScan

## Introduction
`RegionScan` is a R package for comprehensive and scalable region-level genome-wide association testing of alternative region-level multiple-variant and single-variant statistics and visualization of the results. It implements various state-of-the-art region-level tests to improve signal detection under heterogeneous genetic architectures and comparison of multiple-variant region-level and single-variant test results. It leverages LD-based genomic partitioning for LD-adaptive region definition. `RegionScan` is compatible with VCF input file format, and accommodates parallel region-level processing and analysis to improve computational time and memory efficiency. It also provide options for analysis of multi-allelic variants, unbalanced binary phenotypes, with detailed outputs for results interpretation, and utility functions for visualization, comparison, and interpretation.

The detailed information about `RegionScan` can be found in our paper 'link to biorxiv paper??'.

## Download and install
Use the `devtools` package (available from
[CRAN](http://cran-r.c3sl.ufpr.br/web/packages/devtools/index.html)) to
install automatically from this GitHub repository:

```{r, eval=TRUE}
library(devtools)
install_github("brossardMyriam/RegionScan")
```

## Usage
Please see the *add vignette*

## Overview of `RegionScan` capabilities
Overview of RegionScan for genome-wide region-level association analysis on chr19 (to facilitate visualization) in a Usage case study (see details in our [preprint](https://www.biorxiv.org/content/10.1101/2024.03.04.582374v1)). 
![image](https://github.com/brossardMyriam/RegionScan/assets/49451601/615abbbe-0296-4f4d-8d17-137c09511a98)

## Authors
- Myriam Brossard (brossard@lunenfeld.ca), author and main developer
- Yun J Yoo, contributor

## Citation
RegionScan : A comprehensive R package for region-level genome-wide association testing with integration and visualization of multiple-variant and single-variant hypothesis testing.
Brossard M, Roshandel D, Luo K, Yavartanoo F, Paterson AD, Yoo YJ, Bull SB.
DOI: https://www.biorxiv.org/content/10.1101/2024.03.04.582374v1

## License
This package is released under the [GNU General Public License (GPL) v3.0.](https://www.gnu.org/licenses/gpl-3.0.html)
