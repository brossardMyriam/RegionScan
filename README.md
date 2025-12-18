# RegionScan : A comprehensive R package for region-level genome-wide association testing with integration and visualization of multiple-variant and single-variant hypothesis testing 

## Introduction
`RegionScan` is a R package for comprehensive and scalable region-level genome-wide association testing of alternative region-level multiple-variant and single-variant statistics and visualization of the results. It implements various state-of-the-art region-level tests to improve signal detection under heterogeneous genetic architectures and comparison of multiple-variant region-level and single-variant test results. It leverages LD-based genomic partitioning for LD-adaptive region definition. `RegionScan` is compatible with VCF input file format, and accommodates parallel region-level processing and analysis to improve computational time and memory efficiency. It also provide options for analysis of multi-allelic variants, unbalanced binary phenotypes, with detailed outputs for results interpretation, and utility functions for visualization, comparison, and interpretation.

The detailed information about `RegionScan` can be found in our [paper](https://doi.org/10.1093/bioadv/vbaf052). 


## Download and install
Use the `devtools` package (available from
[CRAN](http://cran-r.c3sl.ufpr.br/web/packages/devtools/index.html)) to
install automatically from this repository:

```{r, eval=TRUE}
library(devtools)
install_github("brossardMyriam/RegionScan", build_vignettes = TRUE)
```

## Usage
See the [tutorial](https://github.com/brossardMyriam/RegionScan/blob/main/vignettes/RegionScan.pdf) for an illustration on how to use `RegionScan` on the realistic example dataset provided with the package.

## Overview of `RegionScan` capabilities
Overview of RegionScan for genome-wide region-level association analysis on chr19 (to facilitate visualization) in a Usage case study (see details in [Bioinformatics Advances paper](https://doi.org/10.1093/bioadv/vbaf052)). 
![image](https://github.com/brossardMyriam/RegionScan/assets/49451601/615abbbe-0296-4f4d-8d17-137c09511a98)

## Contributors
- Myriam Brossard (brossard@lunenfeld.ca), author and main developer
- Yun J Yoo, contributor

## Citation
RegionScan: A comprehensive R package for region-level genome-wide association testing with integration and visualization of multiple-variant and single-variant hypothesis testing.
Brossard M, Roshandel D, Luo K, Yavartanoo F, Paterson AD, Yoo YJ, Bull SB.
Bioinformatics Advances, vbaf052, https://doi.org/10.1093/bioadv/vbaf052

## License
This package is released under the [GNU General Public License (GPL) v3.0.](https://www.gnu.org/licenses/gpl-3.0.html)

## Dashboard - NEW !!!
An online dashboard for **exploring and visualizing** RegionScan outputs is now available. Simply upload your RegionScan output files, and the dashboard will generate **interactive visualizations**, including:
QQ plots, Miami plots, locus-level plots with options for gene annotations and recombination rate displays—as well, as well as tools for filtering results and exporting tables. These features are designed to facilitate investigation and comparison of results across methods, regions, and more! The dashboard is implemented in a flexible way and can also accommodate additional region-level tests beyond those currently implemented in RegionScan, provided the input files follow a similar format. For more details, please visit the [RegionScan dashboard page!](https://bulllab-lunenfeld.github.io/region-scan-dashboard/)
![image](https://github.com/brossardMyriam/RegionScan/blob/main/RegionScan_dashboard.png)
