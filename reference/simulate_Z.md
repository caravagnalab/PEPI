# Synthetic generative model for the unified multirates model

Simulates a deterministic mean-field trajectory (the same 2x2
growth/switch ODE the Stan model itself uses as its likelihood mean) for
a wt-background lineage plus any number of driver, driver-only,
driver+clade, and wt-background-clade sub-lineages, each founded by a
single cell at its own introduction time. CCF/fraction/mutation-count
observations are then generated with noise matched to the Stan model's
own likelihoods (beta_proportion for CCF/fractions, lognormal for
counts, Poisson with the same rate formulas as the model's \`generated
quantities\` block for mutation counts), so that data simulated here is
recoverable by \`fit_multirates()\`. Integrate the 2-epistate
growth/switch ODE analytically via matrix exponential.

## Usage

``` r
simulate_Z(lambda_n, s, omega_n, omega_p, t1, t2, Z0)
```

## Arguments

- lambda_n:

  Growth rate in the "-" epistate.

- s:

  Fitness of "+" relative to "-" (lambda_p = lambda_n\*(1+s)).

- omega_n, omega_p:

  Switch rates - to + and + to -.

- t1, t2:

  Start/end time.

- Z0:

  Length-2 initial state c(z_minus, z_plus) at t1.

## Value

Length-2 state c(z_minus, z_plus) at t2.
