# Observed data in the format of `get_predicted_trajectories()`.

Counts of each population in the "-"/"+" samples are its CCF times the
size of the sample, as in the Stan model's transformed data; total
counts at the intermediate times are added to \`pop = "total"\`.

## Usage

``` r
.observed_trajectories(x)
```

## Value

A tibble with columns \`time\`, \`pop\`, \`id\`, \`z_neg\`, \`z_pos\`,
\`z_tot\`, \`frac_pos\`, \`ccf_neg\`, \`ccf_pos\`.
