# Build a default set of initial values for the multirates model.

Several time parameters have bounds that depend on the \*sampled\* value
of \`tmrca\` (e.g. \`t_clade_wt\` is bounded above by the earliest
sampling time, below by \`tmrca\`), so Stan's default random
initialization frequently draws an infeasible \`tmrca\` and fails before
sampling starts. This builds a single feasible starting point instead.

## Usage

``` r
.multirates_default_init(data)
```

## Arguments

- data:

  Stan data list, as built by
  [`build_stan_data_multirates()`](https://caravagnalab.github.io/PEPI/reference/build_stan_data_multirates.md).

## Value

A list (of length 1) suitable for cmdstanr's \`init\` argument.
