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
sim = simulate_multirates_tree(
  sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
  lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
  clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
  mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
)
x = init_multirates(sim$tables, m_trunk = sim$truth$m_trunk)
```
