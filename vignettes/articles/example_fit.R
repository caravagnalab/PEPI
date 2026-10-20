# Shared by Multirates_inference.Rmd and Convergence_diagnostics.Rmd
# (read with knitr::read_chunk), so that both articles use the same fit.

## ---- prepare-data ----

data("counts")
data("muts")
data("forest")
data("mu")

# The two sampling times.
sampling_times = c(
  counts$sampling_time_1 %>% unique(),
  counts$sampling_time_2 %>% unique()
)
t1 = sampling_times[1]
t2 = sampling_times[2]

# Fraction of "-" and "+" cells in the two samples.
ccf = forest %>% filter(!is.na(sample)) %>% group_by(epistate, sample) %>%
  summarize(n = length(cell_id)) %>% group_by(sample) %>% mutate(n = n/sum(n))

# Total population size at the two sampling times...
z1 = counts %>% filter(time == t1) %>%
  group_by(epistate, mutant) %>% summarize(count = max(count)) %>%
  pull(count) %>% unique() %>% sum()

z2 = counts %>% filter(time == t2) %>%
  group_by(epistate, mutant) %>% summarize(count = max(count)) %>%
  pull(count) %>% unique() %>% sum()

set.seed(31)

# ... measured with lognormal noise.
z1 = exp(rnorm(n = 1, mean = log(z1), sd = 0.2))
z2 = exp(rnorm(n = 1, mean = log(z2), sd = 0.2))

zminus = c(z1, z2) * (ccf %>% filter(epistate == "-") %>% pull(n))
zplus = c(z1, z2) * (ccf %>% filter(epistate == "+") %>% pull(n))

# The driver is the first cell of mutant "B": its branch carries the
# mutations accumulated from the root to the driver acquisition.
driver_id = forest %>% filter(mutant == "B") %>% pull(cell_id) %>% min()
driver_branch = PEPI:::get_branch(forest, cell_id = driver_id, stop_id = 0)

m_driver = muts %>% filter(first_cell_id %in% driver_branch) %>% nrow()
m_driver = rnbinom(n = 1, mu = m_driver, size = 100)

# Driver CCF in the four sorted samples, from the VAF of the driver
# mutations (VAF 1-, VAF 1+, VAF 2-, VAF 2+).
ccf_dr = 2 * (muts %>% filter(first_cell_id %in% driver_id) %>%
  dplyr::select(Sampling_1_N.VAF, Sampling_1_P.VAF,
                Sampling_2_N.VAF, Sampling_2_P.VAF) %>%
  reshape2::melt(id.vars = NULL) %>%
  group_by(variable) %>% summarize(v = mean(value)) %>% pull(v))

## ---- driver-prior ----

# Gamma prior on the "+" -> "-" switch rate of the driver, from a plausible
# range of values. The growth rate is bounded above by the growth of the
# total population (log-linear fit of noisy counts above 1e5 cells) and below
# by the growth of the driver "+" cells between the two samples.
z = counts %>% group_by(time) %>% summarize(tot = sum(count)) %>%
  filter(tot > 1e5)

z$tot = exp(rnorm(n = length(z$tot), mean = log(z$tot), sd = 0.2))

fit = lm(log(tot) ~ time, data = z)

lambda_max = coef(fit)[["time"]]
lambda_min = log((zplus[2] * ccf_dr[4]) / (zplus[1] * ccf_dr[2])) / ((1 + 2) * (t2 - t1))

s_max = 2
s_min = 0.1

# ratio of driver "-" and "+" cells at the second sample
rho = (zminus[2] * ccf_dr[3]) / (zplus[2] * ccf_dr[4])

omega_max = rho * lambda_min * s_max
omega_min = rho * lambda_max * s_min

omega_gamma = PEPI:::fit_gamma_quantiles(omega_min, omega_max, p1 = 0.05, p2 = 0.95)

## ---- tables ----

tables = list(

  population_sizes = tibble(
    time = sampling_times,
    type = "sampling",
    zminus = zminus,
    zplus = zplus
  ),

  # the wild type is the complement of the driver
  wt_ccf = tibble(
    time = rep(sampling_times, each = 2),
    epistate = rep(c("-", "+"), 2),
    ccf = 1 - ccf_dr
  ),

  driver_muts = tibble(
    driver_id = "B",
    m_driver = m_driver,
    ms_driver = -0.8,
    sigma_driver = 0.9,
    alpha_n_driver = omega_gamma$shape,
    beta_n_driver = omega_gamma$rate,
    alpha_p_driver = 0.3,
    beta_p_driver = 6
  ),

  driver_ccf = tibble(
    driver_id = "B",
    time = rep(sampling_times, each = 2),
    epistate = rep(c("-", "+"), 2),
    ccf = ccf_dr
  )
)

x = init_multirates(tables, m_trunk = 0)

## ---- fit ----

x = fit_multirates(
  x,
  seed = 31,
  chains = 4,
  parallel_chains = 4,
  ndraws = 5000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  init = rep(list(list(omega_p_wt = 3e-3, omega_p_driver = 3e-3)), 4),
  mu = mu$mean,
  l = 51e06,
  t_min = -20,
  ms_epi = -0.8,
  sigma_epi = 0.9,
  alpha_lambda = 10,
  beta_lambda = 10,
  alpha_n_wt = 1,
  beta_n_wt = 1e4,
  alpha_p_wt = 0.2,
  beta_p_wt = 10,
  min_kappa = 3000,
  max_kappa = 5000,
  min_sigma_count = 0.05,
  max_sigma_count = 0.15,
  include_poisson = 1,
  include_trunk = 0
)
