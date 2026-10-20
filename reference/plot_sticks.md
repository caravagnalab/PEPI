# Plot a phylogenetic forest with nodes labelled by mutation class.

Nodes carrying a label (e.g. "driver", "truncal") are highlighted; all
other nodes are labelled "subclonal".

## Usage

``` r
plot_sticks(forest, labels, cls = NULL)
```

## Arguments

- forest:

  Table of forest nodes with columns \`cell_id\`, \`ancestor\`,
  \`mutant\`, \`epistate\`, \`sample\` and \`birth_time\`.

- labels:

  Table with columns \`cell_id\` and \`label\`.

- cls:

  Named vector of label colors. If NULL a palette is generated.

## Value

A ggplot object.

## Examples

``` r
data("forest")
set.seed(31)
cells = sample(forest$cell_id[!is.na(forest$sample)], 20)
ids = unique(c(0, unlist(lapply(cells, function(id) PEPI:::get_branch(forest, id, 0)))))
driver_id = min(forest$cell_id[forest$mutant == "B"])
labels = data.frame(cell_id = PEPI:::get_branch(forest, driver_id, 0), label = "driver")
PEPI:::plot_sticks(forest[forest$cell_id %in% ids, ], labels,
  cls = c(subclonal = "gainsboro", driver = "forestgreen"))
#> Warning: Removed 1 row containing non-finite outside the scale range
#> (`stat_edge_link()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
```
