# Plot a phylogenetic forest coloured by species.

Cells are placed by birth time, coloured by species (mutant + epistate)
and shaped by sample.

## Usage

``` r
plot_forest(
  forest,
  ids = NULL,
  colors = c(`wt-` = "#82CCE0", `wt+` = "#CE6D6E", `B-` = "steelblue", `B+` =
    "darkorange3")
)
```

## Arguments

- forest:

  Table of forest nodes with columns \`cell_id\`, \`ancestor\`,
  \`mutant\`, \`epistate\`, \`sample\` and \`birth_time\`.

- ids:

  Optional vector of cell ids to plot (e.g. the branches of a subset of
  sampled cells). If NULL, the whole forest is plotted.

- colors:

  Named vector of species colors.

## Value

A ggplot object.

## Examples

``` r
data("forest")
set.seed(31)
cells = sample(forest$cell_id[!is.na(forest$sample)], 20)
ids = unique(c(0, unlist(lapply(cells, function(id) PEPI:::get_branch(forest, id, 0)))))
PEPI:::plot_forest(forest, ids)
```
