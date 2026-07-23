# Associate colors to nodes.

A list of colors labelled by nodes is generated.

## Usage

``` r
get_colors(max_depth)
```

## Arguments

- max_depth:

  Maximal number of levels of the tree

## Value

Named list of colors.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
PEPI:::get_colors(max_depth = 2)
#>           -          --          -+         ---         --+         -+- 
#> "#5050FFFF" "#CE3D32FF" "#749B58FF" "#F0E685FF" "#466983FF" "#BA6338FF" 
#>         -++ 
#> "#5DB1DDFF" 
```
