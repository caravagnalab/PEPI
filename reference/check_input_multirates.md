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
if (FALSE) { # \dontrun{
check_input_multirates(tables)
} # }
```
