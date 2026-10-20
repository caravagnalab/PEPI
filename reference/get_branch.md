# Get the branch of a cell in a phylogenetic forest.

Walks the forest from \`cell_id\` back through its ancestors, stopping
before \`stop_id\`.

## Usage

``` r
get_branch(forest, cell_id, stop_id)
```

## Arguments

- forest:

  Table of forest nodes with columns \`cell_id\` and \`ancestor\` (as
  returned by ProCESS' `forest$get_nodes()`).

- cell_id:

  Id of the cell the branch starts from.

- stop_id:

  Id of the ancestor where the walk stops (excluded).

## Value

Vector of cell ids, from \`cell_id\` up to the child of \`stop_id\`.

## Examples

``` r
data("forest")
driver_id = min(forest$cell_id[forest$mutant == "B"])
PEPI:::get_branch(forest, cell_id = driver_id, stop_id = 0)
#> [1] 57 36 27 13  8  3  1
```
