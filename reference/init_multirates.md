# Create a PEPI object for the unified multirates model.

A PEPI object holding the tables required to build the data list for
`inst/multirates_positive_s.stan` is created. See
[`check_input_multirates()`](https://caravagnalab.github.io/PEPI/reference/check_input_multirates.md)
for the required table schema.

## Usage

``` r
init_multirates(tables, m_trunk = 0L)
```

## Arguments

- tables:

  Named list of tibbles: `population_sizes`, `wt_ccf` (required), and
  any subset of `driver_muts`/`driver_ccf`,
  `driver_n_muts`/`driver_n_ccf`, `driver_p_muts`/`driver_p_ccf`,
  `dc_muts`/`dc_ccf`, `clades_wt`/`clades_wt_ccf` (optional).

- m_trunk:

  Number of truncal mutations.

## Value

PEPI object of class "PEPI_Multirates"

## Examples

``` r
if (FALSE) { # \dontrun{
init_multirates(tables, m_trunk = 120)
} # }
```
