# Plot the simulated data: phylogenetic forest and VAF spectra.

Samples \`n_cells\` cells from each sample and restricts the forest to
their branches. Top: the forest coloured by species
([`plot_forest()`](https://caravagnalab.github.io/PEPI/reference/plot_forest.md))
and with the driver branch and the truncal node highlighted
([`plot_sticks()`](https://caravagnalab.github.io/PEPI/reference/plot_sticks.md)).
Bottom: pairwise and marginal VAF spectra of the sorted samples
([`plot_vaf()`](https://caravagnalab.github.io/PEPI/reference/plot_vaf.md)),
with mutations labelled as truncal, driver (acquired on the branch from
the root to the first cell of \`driver_mutant\`) or subclonal.

## Usage

``` r
plot_data(
  forest,
  muts,
  n_cells = 50,
  seed = NULL,
  driver_mutant = "B",
  vaf_cut = 0.03,
  cls = c(subclonal = "gainsboro", driver = "forestgreen", truncal = "darkgoldenrod3")
)
```

## Arguments

- forest:

  Table of forest nodes (see
  [`plot_forest()`](https://caravagnalab.github.io/PEPI/reference/plot_forest.md)).

- muts:

  Mutation table with VAF columns and \`first_cell_id\`.

- n_cells:

  Number of cells drawn from each sample.

- seed:

  Optional seed for the cell sampling.

- driver_mutant:

  Name of the driver mutant.

- vaf_cut:

  Minimum VAF shown in the marginal histograms.

- cls:

  Named vector of colors for the labels \`subclonal\`, \`driver\`,
  \`truncal\`.

## Value

A ggplot object.

## Examples

``` r
data("forest")
data("muts")
PEPI:::plot_data(forest, muts, seed = 31)
#> Warning: Removed 1 row containing non-finite outside the scale range
#> (`stat_edge_link()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
```
