# Plot marginal and pairwise VAF distributions of the sorted samples.

Every column whose name contains "VAF" is treated as one sample (columns
starting with "normal" are dropped). Columns named as returned by
ProCESS' `simulate_seq()` with an epistate labelling (e.g.
\`Sampling_1_N.VAF\`, \`Sampling_1_P.VAF\`) are renamed to \`VAF 1-\`,
\`VAF 1+\`, ...

## Usage

``` r
plot_vaf(muts, colors = NULL, vaf_cut = 0.02, size = 0.01, assemble = TRUE)
```

## Arguments

- muts:

  Mutation table with VAF columns and, optionally, a \`label\` column.

- colors:

  Named vector of label colors.

- vaf_cut:

  Minimum VAF shown in the marginal histograms.

- size:

  Point size of the pairwise scatter plots.

- assemble:

  If TRUE, return one figure with the pairwise scatter plots on top and
  the marginal histograms below; otherwise a list with elements
  \`marginals\` and \`pairs\`.

## Value

A ggplot object, or a list of ggplot objects.

## Examples

``` r
data("muts")
PEPI:::plot_vaf(muts, vaf_cut = 0.03, size = 0.1)
```
