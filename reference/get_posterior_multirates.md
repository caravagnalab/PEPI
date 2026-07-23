# Extract and tidy the full posterior (and built-in prior draws) from a multirates fit.

Every recognized Stan variable (raw parameters and all \`generated
quantities\` families: z/frac/ccf/m posterior predictives, plus the
inline \`\*\_prior\` draws) is parsed into one long tidy table with
columns \`type\` (prior/posterior), \`family\` (z/frac/ccf/m/param),
\`group\` (wt/driver/driver_n/driver_p/dc/clade_wt/global), \`id\`,
\`sampling_time\`, \`epistate\`, \`event\`, \`quantity\` (driver/clade
leg for dc entities), \`base\`, \`variable\`, draw indices, and
\`value\`.

## Usage

``` r
get_posterior_multirates(x)
```

## Arguments

- x:

  PEPI_Multirates object with \`x\$inference\$multirates\` populated by
  [`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md).

## Value

\`x\` with \`x\$posterior\$multirates\` populated.

## Examples

``` r
if (FALSE) { # \dontrun{
get_posterior_multirates(x)
} # }
```
