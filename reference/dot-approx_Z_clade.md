# R version of the Stan approx_Z_clade(), vectorised over posterior draws.

Clades start from a single "+" cell (boundary condition), so there is no
chi term.

## Usage

``` r
.approx_Z_clade(lambda_n, lambda_p, omega_n, omega_p, t_1, t_2, W1)
```

## Value

A list with the "-" (\`neg\`) and "+" (\`pos\`) clade sizes.
