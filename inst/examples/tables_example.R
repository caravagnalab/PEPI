# Input tables of the polyclonal example shipped with PEPI (simulated with
# ProCESS, see the "Data generation" article), as in the "Inference" article.

library(PEPI)
library(dplyr)

data("counts", package = "PEPI")
data("muts", package = "PEPI")
data("forest", package = "PEPI")
data("mu", package = "PEPI")

sampling_times = c(unique(counts$sampling_time_1), unique(counts$sampling_time_2))

ccf = forest %>% filter(!is.na(sample)) %>% count(epistate, sample) %>%
  group_by(sample) %>% mutate(n = n/sum(n)) %>% ungroup()

z = sapply(sampling_times, function(t){
  counts %>% filter(time == t) %>% group_by(epistate, mutant) %>%
    summarize(count = max(count), .groups = "drop") %>% pull(count) %>% unique() %>% sum()
})

set.seed(31)
z = exp(rnorm(n = 2, mean = log(z), sd = 0.2))

zminus = z * (ccf %>% filter(epistate == "-") %>% pull(n))
zplus = z * (ccf %>% filter(epistate == "+") %>% pull(n))

driver_id = forest %>% filter(mutant == "B") %>% pull(cell_id) %>% min()
driver_branch = PEPI:::get_branch(forest, cell_id = driver_id, stop_id = 0)

m_driver = rnbinom(n = 1, mu = sum(muts$first_cell_id %in% driver_branch), size = 100)

# driver CCF in the samples 1-, 1+, 2-, 2+
ccf_dr = 2 * colMeans(muts[muts$first_cell_id == driver_id,
                           c("Sampling_1_N.VAF", "Sampling_1_P.VAF", "Sampling_2_N.VAF", "Sampling_2_P.VAF")])

tables = list(
  population_sizes = tibble(time = sampling_times, type = "sampling", zminus = zminus, zplus = zplus),
  wt_ccf = tibble(time = rep(sampling_times, each = 2), epistate = rep(c("-", "+"), 2), ccf = 1 - ccf_dr),
  driver_muts = tibble(driver_id = "B", m_driver = m_driver, ms_driver = -0.8, sigma_driver = 0.9,
                       alpha_n_driver = 1.5, beta_n_driver = 4, alpha_p_driver = 0.3, beta_p_driver = 6),
  driver_ccf = tibble(driver_id = "B", time = rep(sampling_times, each = 2),
                      epistate = rep(c("-", "+"), 2), ccf = ccf_dr)
)

x = init_multirates(tables, m_trunk = 0)
