# PEPI

The `PEPI` package implements a Bayesian statistical framework to
quantify epigenetic plasticity from genomics data. The model considers
an experimental design in which cells are sorted according to their
epimutation state (positive versus negative) prior to sequencing;
Variant Allele Frequencies (VAFs) are then obtained from the two sorted
populations, together with the corresponding cell counts, across driver
and wild-type clades. This procedure can be repeated across multiple
time points, enabling the study of epigenetic and clonal dynamics over
time via a unified branching-process model fit with Stan.

#### Citation

A manuscript describing `PEPI` is in preparation; there is no preprint
yet. In the meantime, if you use `PEPI` please cite this repository.

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/PEPI/-yellow.svg)](https://caravagnalab.github.io/PEPI)

## Installation

You can install the released version of **PEPI** from
[GitHub](https://github.com) with:

``` r

devtools::install_github("caravagnalab/PEPI")
```

**Note:** PEPI requires [cmdstanr](https://mc-stan.org/cmdstanr/) to fit
models. Please follow the [installation
instructions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
before use.

------------------------------------------------------------------------

#### Copyright and contacts

Cancer Data Science (CDS) Laboratory. University of Trieste, Italy.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
