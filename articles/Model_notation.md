# Model and notation

This page works through the mathematical model implemented by
`inst/multirates_positive_s.stan`, and lays each quantity next to the R
argument, table or column that supplies it. Read this before
`vignette("Multirates_inference")` if you want to know *why* the input
tables are shaped the way they are, rather than just how to build them.

## Two epigenetic states, growth and switching

Cells exist in one of two epigenetic states, written `-` and `+`. A
population’s size in each state is a vector $`Z(t) = (Z_-(t), Z_+(t))`$,
and every lineage in the model (the wild-type background, a driver
clade, a driver+clade combo, …) is governed by the same linear ODE, just
with its own rate parameters:

``` math
\frac{dZ}{dt} = A\,Z, \qquad
A = \begin{pmatrix} \lambda_- & \omega_- \\ \omega_+ & \lambda_+ \end{pmatrix}
```

- $`\lambda_-`$ is the division rate of `-` cells.
- $`\lambda_+ = \lambda_-(1+s_+)`$ is the division rate of `+` cells,
  written relative to $`\lambda_-`$ via a fitness coefficient $`s_+`$
  (so $`s_+>0`$ means `+` cells divide faster).
- $`\omega_-`$ and $`\omega_+`$ are the rates at which cells switch
  epigenetic state.

Because $`A`$ does not depend on $`t`$, the population at any later time
is an exact matrix-exponential solution:

``` math
Z(t_2) = \exp\!\big(A\,(t_2-t_1)\big)\,Z(t_1)
```

This is the Stan model’s own `Z()` function.

## The groups: wt, drivers, driver-only clades, driver+clade combos

The population is not a single lineage: it is the wild-type background
plus any number of sub-lineages, each founded at its own time from a
single cell and governed by the ODE above with its own rate parameters.

| Genotype | Founded at | `s_g` | Switch rates | Confined to one state? |
|----|----|----|----|----|
| **wt** (background) | $`t_{mrca}`$ | — | $`\omega_-, \omega_+`$ | no |
| **clade_wt** (wt sub-lineage) | $`t_{clade}\geq t_{mrca}`$ | — | $`\omega_-, \omega_+`$ (shared with wt) | no |
| **driver** | $`t_{driver}\geq t_{mrca}`$ | $`s_{driver_i}`$ (own) | $`\omega_-^{driver_i}, \omega_+^{driver_i}`$ (own) | no |
| **driver_n** | $`t_{driver^{-}}\geq t_{mrca}`$ | $`s_{driver^{-}_i}`$ (own) | none (no switching) | yes, `-` only |
| **driver_p** | $`t_{driver^{+}}\geq t_{mrca}`$ | $`s_{driver^{+}_i}`$ (own) | none (no switching) | yes, `+` only |
| **dc** driver leg | $`t_{dc,1} = t_{mrca}+\delta_{t,1}`$ | $`s_{dc_i}`$ (own) | $`\omega_-^{dc_i}, \omega_+^{dc_i}`$ (own) | no |
| **dc** clade leg | $`t_{dc,2} = t_{dc,1}+\delta_{t,2}`$ | $`s_{dc_i}`$ (shared with its driver leg) | $`\omega_-^{dc_i}, \omega_+^{dc_i}`$ (shared with its driver leg) | no |

whereas `s_g` is genotype specific fitness advantage with respect to the
wild type. `driver_n`,`driver_p` are the simplest cases: with no
switching, the ODE collapses to plain exponential growth,
$`Z(t) = \exp\big(\lambda_-(1+s_{driver^{\pm}})(t-t_{driver^{\pm}})\big)`$.
A `dc` (driver+clade) entity is a *nested* pair: its driver leg is
founded first with a single cell, and its clade leg is a further
sub-lineage founded later *within* the driver leg, sharing the same
growth/switch rates - modelling a clade defined by a further epimutation
switch event that occurred after the driver mutation.

At each sampling time, the total population in each state is the sum of
the wt background and every additive group (`driver`, `driver_n`,
`driver_p`, and the `dc` driver leg - `clade_wt` and the `dc` clade leg
are *nested inside* their parent lineage’s mass, not counted again):

``` math
Z^{tot}_-(t) = Z_{wt,-}(t)+\sum_i Z_{driver_i,-}(t)+\sum_i Z_{driver^{-}_i}(t)+\sum_i Z_{dc_i,-}(t)
```

and analogously for $`Z^{tot}_+(t)`$.

## The observation model

### Cancer cell fractions (CCF) and mutation-arm fractions

For genotypes that occupy both epigenetic states, the CCF of the `-`
phenotype at a sampling time is the fraction of that genotype’s cells
that are in the `-` state:

``` math
CCF_{g,-}(t)=\frac{Z_{g,-}(t)}{Z_{g,-}(t)+Z_{g,+}(t)}
```

This is observed (as data) or predicted (as a posterior draw) with Beta
noise of concentration $`\kappa`$:

``` math
CCF_{g,-}^{obs} \sim \text{Beta}\big(CCF_{g,-}(t)\cdot\kappa,\ (1-CCF_{g,-}(t))\cdot\kappa\big)
```

This CCF is defined only for genotypes `g` with both `-` and `+`
compartments such as `wt` and `driver`; genotypes confined to a single
epigenetic state, such as `driver_n`, `driver_p` do not have a
phenotypic CCF split.

`kappa` is the same idea as `bayesplot`’s `neff_ratio` - it controls how
tightly the observed fraction is expected to hug its predicted mean;
high $`\kappa`$ means a very peaked, low-noise likelihood (see
`vignette("Convergence_diagnostics")` for why a very high `kappa` can
make Stan’s sampler harder to initialize).

### Population sizes

Total population sizes for each genotype at sampling times (`zminus`,
`zplus`) and at any additional intermediate times (`ztot`) are observed
with lognormal noise of scale $`\sigma_{count}`$:

``` math
z_g^{obs}(t) \sim \text{LogNormal}\big(\log z_g(t),\ \sigma_{count}\big)
```

### Mutation counts

Every branch event (the trunk, a driver acquisition, a clade split, …)
accumulates mutations as a Poisson process while cells divide. For a
branch with growth rate $`\lambda`$ (its own $`\lambda_-`$ or
$`\lambda_-(1+s)`$ once a fitness effect has taken hold) spanning an
elapsed time $`\Delta t`$:

``` math
m \sim \text{Poisson}\big(4\,\mu\,\ell\,\lambda\,\Delta t\big)
```

where $`\mu`$ (`mu`) is the per-division, per-bp mutation rate and
$`\ell`$ (`l`) is the genome length. The trunk uses
$`\Delta t = t_{mrca} - t_{min}`$; a driver or `driver_n` acquisition
uses $`\Delta t = t_{driver}-t_{mrca}`$; a `dc` combo’s clade leg uses
the *post-driver* rate $`\lambda_-(1+s_{dc})`$ over
$`\Delta t = t_{dc,2}-t_{dc,1}`$.

## Priors

| Parameter | Prior | Set via |
|----|----|----|
| $`\lambda_-`$ (`lambda_n`) | $`\text{Gamma}(\alpha_\lambda,\beta_\lambda)`$ | `alpha_lambda`, `beta_lambda` |
| $`s_+`$ (`s_epi`) | $`\text{LogNormal}(ms_{epi},\sigma_{epi})`$ | `ms_epi`, `sigma_epi` |
| $`\omega_-,\omega_+`$ | $`\text{Gamma}(\alpha,\beta)`$ | `alpha_n_wt`/`beta_n_wt`, `alpha_p_wt`/`beta_p_wt` |
| $`t_{mrca}`$ | $`\text{Uniform}(t_{min}, t_{max})`$ | `t_min` (`t_max` is the latest sampling time) |
| $`s_{driver_i}, s_{driver^{-}_i}, s_{driver^{+}_i}, s_{dc_i}`$ | $`\text{LogNormal}(ms,\sigma)`$, one prior per entity | `ms_driver`/`sigma_driver` columns (etc.) in the entity’s `_muts` table |
| $`\omega_-^{driver_i},\omega_+^{driver_i}`$ (and `_dc`) | $`\text{Gamma}(\alpha,\beta)`$, one prior per entity | `alpha_n_driver`/`beta_n_driver` columns (etc.) |
| $`\kappa`$ | $`\text{Uniform}(min_\kappa, max_\kappa)`$ | `min_kappa`, `max_kappa` |
| $`\sigma_{count}`$ | $`\text{Uniform}(min_{\sigma},max_{\sigma})`$ | `min_sigma_count`, `max_sigma_count` |

Global priors are arguments to
[`fit_multirates()`](https://caravagnalab.github.io/PEPI/reference/fit_multirates.md).
Per-entity priors (driver/driver_n/driver_p/dc fitness and switch-rate
hyperparameters) are supplied per row in the corresponding `_muts`
table, since a well-studied oncogene and a novel candidate driver may
legitimately warrant different priors within the same fit.

## Notation table

The table below is the authoritative map from math to code: the Stan
parameter/generated-quantity name (as it appears in
[`get_posterior_multirates()`](https://caravagnalab.github.io/PEPI/reference/get_posterior_multirates.md)’s
tidy output), and the R-facing input that supplies it.

| Symbol | Stan / getter name | R input |
|:---|:---|:---|
| $`\lambda_-`$ | `lambda_n` | `fit_multirates(alpha_lambda=, beta_lambda=)` |
| $`s_{+}`$ ($`\lambda_+=\lambda_-(1+s_{+})`$) | `s_epi` | `fit_multirates(ms_epi=, sigma_epi=)` |
| $`\omega_{+}`$ | `omega_p_wt` | `fit_multirates(alpha_p_wt=, beta_p_wt=)` |
| $`\omega_{-}`$ | `omega_n_wt` | `fit_multirates(alpha_n_wt=, beta_n_wt=)` |
| $`t_{mrca}`$ | `tmrca` | inferred (bounded by `t_min` and the latest sampling time) |
| $`t_{min}`$ | `t_min` | `fit_multirates(t_min=)` |
| $`s_{driver_i}`$ | `s_driver` | `driver_muts$ms_driver`, `sigma_driver` |
| $`\omega_{+}^{driver_i}`$, $`\omega_{-}^{driver_i}`$ | `omega_p_driver`, `omega_n_driver` | `driver_muts$alpha_n_driver`, … (per row) |
| $`t_{driver_i}`$ | `t_driver` | inferred (bounded by `tmrca` and `t_max`) |
| $`m_{driver_i}`$ (observed) | `m_driver` | `driver_muts$m_driver` |
| $`s_{driver^{-}_i}`$ / $`s_{driver^{+}_i}`$ | `s_driver_n` / `s_driver_p` | `driver_n_muts`/`driver_p_muts$ms_driver_n`/`_p` |
| $`t_{driver^{-}_i}`$ / $`t_{driver^{+}_i}`$ | `t_driver_n` / `t_driver_p` | inferred |
| $`s_{dc_i}`$ | `s_dc` | `dc_muts$ms_dc`, `sigma_dc` |
| $`t_{dc_i,1}, t_{dc_i,2}`$ | `delta_t_dc[,1]`, `delta_t_dc[,2]` | inferred |
| $`t_{clade_i}`$ (wt sub-lineage) | `t_clade_wt` | inferred |
| $`CCF_{g,-}(t)`$ (observed) | `ccf` (family `"ccf"`) | `driver_ccf`/`clades_wt_ccf`/`dc_ccf`/`wt_ccf$ccf` |
| $`\kappa`$ | `kappa` | `fit_multirates(min_kappa=, max_kappa=)` |
| $`z_-(t), z_+(t)`$ (observed) | `zminus`, `zplus` | `population_sizes$zminus`, `zplus` |
| $`\sigma_{count}`$ | `sigma_count` | `fit_multirates(min_sigma_count=, max_sigma_count=)` |
| $`\mu`$ (mutation rate) | `mu` | `fit_multirates(mu=)` |
| $`\ell`$ (genome length) | `l` | `fit_multirates(l=)` |

For the full input table schemas (columns required per group), see
[`?init_multirates`](https://caravagnalab.github.io/PEPI/reference/init_multirates.md)
and
[`?check_input_multirates`](https://caravagnalab.github.io/PEPI/reference/check_input_multirates.md),
or the worked example in `vignette("Multirates_inference")`.
