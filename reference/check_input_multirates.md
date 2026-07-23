# Check if multirates input data has the correct form.

Validates a named list of tibbles against the schema required to build
the data list for `inst/multirates_positive_s.stan`: required tables and
columns are present, every id referenced in a \`\_ccf\` table exists in
its metadata table (and vice versa), every \`time\` referenced in a
\`\_ccf\`/ \`wt_ccf\` table exists in \`population_sizes\` with the
right \`type\`, and every id x time (x epistate) combination is present
exactly once.

## Usage

``` r
check_input_multirates(tables)
```

## Arguments

- tables:

  Named list of tibbles, see
  [`init_multirates()`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md).

## Value

Invisibly TRUE if valid, otherwise stops with a descriptive message.

## Examples

``` r
sim = simulate_multirates_tree(
  sampling_times = c(3, 6, 9), tmrca = 1, t_min = 0,
  lambda_n = 1, s_epi = 0.15, omega_n_wt = 5e-3, omega_p_wt = 2e-3,
  clades_wt = list(list(id = "c1", t_clade_wt = 1.5)),
  mu = 1e-7, l = 2.7e9, kappa = 20, sigma_count = 0.1, seed = 42
)
check_input_multirates(sim$tables)
```
