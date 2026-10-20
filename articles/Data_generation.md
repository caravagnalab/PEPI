# Data generation

The datasets shipped with PEPI (`counts`, `forest`, `muts` and `mu`) are
a polyclonal tumour simulated with
[ProCESS](https://caravagnalab.github.io/ProCESS/)
([GitHub](https://github.com/caravagnalab/ProCESS)): a wild-type
population that switches between the “-” and “+” epistates, in which a
driver (mutant `B`) is acquired by a “-” cell. This article shows how
they were generated. The code is not run here: simulating the growth up
to 3 million cells, sequencing the samples and fitting the mutation rate
take a while.

Install ProCESS with `remotes::install_github("caravagnalab/ProCESS")`.

``` r

library(ProCESS)
library(dplyr)
```

## Parameters

``` r

lambda_n = 1        # growth rate of "-" cells
s_epi = 0.5         # fitness of "+" cells: lambda_p = lambda_n * (1 + s_epi)
omega_p = 8e-3      # switch rate "-" -> "+"
omega_n = 1e-4      # switch rate "+" -> "-"

t_dr = 4            # time of the driver acquisition
s_dr = 0.35         # fitness of the driver
omega_p_dr = 4e-2   # driver switch rate "-" -> "+"
omega_n_dr = 0.3    # driver switch rate "+" -> "-"

set.seed(63)
```

## Tissue simulation

A single wild-type “-” cell grows without death. At time `t_dr` a random
“-” cell acquires the driver, which multiplies the growth rate of both
epistates by `1 + s_dr` and has its own switch rates.

``` r

sim = SpatialSimulation("polyclonal", width = 1e4, height = 1e4, save_snapshot = TRUE)

sim$border_growth_model = FALSE
sim$history_delta = 0.1

sim$add_mutant(name = "wt",
               epigenetic_rates = c("+-" = omega_n, "-+" = omega_p),
               growth_rates = c("+" = lambda_n * (1 + s_epi), "-" = lambda_n),
               death_rates = c("+" = 0.0, "-" = 0.0))

sim$place_cell("wt-", 5e3, 5e3)

sim$run_up_to_time(t = t_dr)

sim$add_mutant(name = "B",
               epigenetic_rates = c("+-" = omega_n_dr, "-+" = omega_p_dr),
               growth_rates = c("+" = lambda_n * (1 + s_epi) * (1 + s_dr),
                                "-" = lambda_n * (1 + s_dr)),
               death_rates = c("+" = 0.0, "-" = 0.0))

candidates = sim$get_cells() %>% filter(epistate == "-")
cell = candidates[sample(x = 1:nrow(candidates), size = 1), ]

sim$mutate_progeny(cell, "B")
```

The tumour is sampled (3000 cells) when it reaches 1 million and 3
million cells.

``` r

wtp = sim$var("wt+")
wtn = sim$var("wt-")
dcp = sim$var("B+")
dcn = sim$var("B-")

sim$run_until(wtp + wtn + dcp + dcn > 1e6)
t1 = sim$get_clock()
sim$sample_cells("Sampling_1", num_of_cells = 3000)

sim$run_until(wtp + wtn + dcp + dcn > 3e6)
t2 = sim$get_clock()
sim$sample_cells("Sampling_2", num_of_cells = 3000)

forest = sim$get_samples_forest()
```

The cell counts over time and the phylogenetic forest of the sampled
cells:

``` r

counts = sim$get_count_history() %>%
  mutate(sampling_time_1 = t1, sampling_time_2 = t2)

forest_nodes = forest$get_nodes()
```

## Mutations and sequencing

Passenger SNVs are placed on the forest with the same rate in both
epistates and in both mutants.

``` r

m_engine = MutationEngine(setup_code = "demo", context_sampling = 5)

m_engine$add_mutant(mutant_name = "wt",
                    passenger_rates = list("+" = c(SNV = 5e-7, CNA = 0),
                                           "-" = c(SNV = 5e-7, CNA = 0)),
                    driver = c())

m_engine$add_mutant(mutant_name = "B",
                    passenger_rates = list("+" = c(SNV = 5e-7, CNA = 0),
                                           "-" = c(SNV = 5e-7, CNA = 0)),
                    driver = c())

m_engine$add_exposure(c(SBS5 = 1))
m_engine$add_exposure(c(ID2 = 1))

phylo_forest = m_engine$place_mutations(forest, 500, 100)
```

Cells are sorted by epistate before sequencing: labelling each cell by
its epistate makes `simulate_seq()` return, for every mutation, the
reads of the “-” (`N`) and “+” (`P`) cells of every sample, at 100x
coverage.

``` r

epi_labelling = function(cell){
  if (cell$epistate == "+") return("P")
  if (cell$epistate == "-") return("N")
  return("")
}

seq_results = simulate_seq(phylo_forest, coverage = 100, cell_labelling = epi_labelling)
```

We keep somatic SNVs and annotate each one with the first cell where it
occurs, which tells on which branch of the forest it was acquired.

``` r

muts = seq_results$mutations %>%
  filter(causes != "germinal", nchar(ref) == 1 & nchar(alt) == 1) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(first_cell_id = phylo_forest$get_first_occurrences(SNV(chr, chr_pos, alt, ref))[[1]]) %>%
  ungroup()
```

## Mutation rate

The mutation rate per division is estimated from the neutral tail of the
VAF spectrum of the first “-” sample with MOBSTERh, and its posterior is
obtained with `mu_posterior()`. Both functions belong to
[MOBSTERh](https://github.com/caravagnalab/mobster/tree/mobsterh), the
`mobsterh` branch of [mobster](https://caravagnalab.github.io/mobster/),
which can be installed with
`remotes::install_github("caravagnalab/mobster@mobsterh")`.

``` r

library(mobster)

data = muts %>%
  dplyr::rename(NV = Sampling_1_N.occurrences,
                DP = Sampling_1_N.coverage,
                VAF = Sampling_1_N.VAF,
                from = chr_pos) %>%
  filter(VAF > 0) %>%
  mutate(to = from + 1, karyotype = "1:1",
         segment_id = "chr22:1:51e06:1:1:1")

fit = mobsterh_fit(
  data,
  lr = 0.01,
  purity = 0.99,
  subclonal_clusters = 0:1,
  N_MAX = 1e6,
  tail = TRUE,
  NV_filter = 5,
  truncate_pareto = FALSE,
  enforce_QC_PASS = FALSE,
  n_t = 10,
  vaf_filter = 0.05
)

mu = mobster:::mu_posterior(fit$best, ncells = 2, lq = 0.01, uq = 0.99)
```

`counts`, `forest_nodes` (saved as `forest`), `muts` and `mu` are the
datasets used in `vignette("Multirates_inference")`.
