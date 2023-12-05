


# Fit a multivariate vaf spectrum with epigenetic tree model.
#
# A cmdstanr fit is returned.
#
# @param spectrum Dataframe with number of variants and depth for any mutation and sample
# @param path_to_model String specifying the path where we wannt to save the stan model for a given depth
# @param cmdstan_path String specifying the path to cmdstan folder
# @param max_depth Maximum number of levels
# @param ndraws Number of draws from the posterior
# @param init List of initialization parameters
# @param seed Seed of the computation
# @param mu Mutation rate per division per bp per allele
# @param l length of the genome
#@ param rho_n Purity of - sample
#@ param rho_p Purity of + sample
# @param nu_t Mean of the beta prior on fraction of truncal mutations
# @param qt Number of trials of the beta prior on fraction of truncal mutations
# @param rate_n Mean of the beta prior on the epimutation rate from - to +
# @param qn Number of trials for the beta prior on the epimutation rate from - to +
# @param rate_p Mean of the beta prior on the epimutation rate from + to -
# @param qp Number of trials for the beta prior on the epimutation rate from - to +
# @param k Number of trials for the beta prior on cluster centroids
# @param gamma Concentration of a Dirichlet distribution to split mutations at any node
# @return cdmstanr fit
# @examples
# fit_tree = function(spectrum,path_to_model = "models",cmdstan_path = "my_cmdstan/",
# max_depth = 2,ndraws = 1000,init = NULL,seed = 15,
# mu = 1e-7,l = 2.7*10^9,rho_n = 1,rho_p = 1,nu_t = 0.1,
# qt = 1e4,rate_n = 1e-3,qn = 1e4,rate_p = 1e-3,
# qp = 1e4,k = 1e4,gamma = 150)
# @export

fit_tree = function(spectrum,path_to_model,cmdstan_path,
                    max_depth = 2,ndraws = 1000,init = NULL,seed = 15,
                    mu = 1e-7,l = 2.7*10^9,rho_n = 1,rho_p = 1,nu_t = 0.1,
                    qt = 1e4,rate_n = 1e-3,qn = 1e4,rate_p = 1e-3,
                    qp = 1e4,k = 1e4,gamma = 150){
  
  dir.create(path_to_model)
  cmdstanr::set_cmdstan_path(cmdstan_path)
  
  model = tree_inference_code(max_depth = max_depth,likelihood = T)
  
  write_stan_file(
    model,
    dir = ".",
    basename = paste0(path_to_model,"/tree_inference_depth_",max_depth,".stan"),
    force_overwrite = FALSE,
    hash_salt = ""
  )
  
  data = list(
    delta_m_n = spectrum %>% nrow(),
    mu = mu,
    l = l,
    nu_t = nu_t,
    qt = qt,
    rate_n = rate_n,
    qn = qn,
    rate_p = rate_p,
    qp = qp,
    n = nrow(spectrum),
    Nn = spectrum$Nx,
    Np = spectrum$Ny,
    DPn = spectrum$DPx,
    DPp = spectrum$DPy,
    gamma = gamma,
    rho_n = rho_n,
    rho_p = rho_p,
    k = k)
  
  file = paste0(path_to_model,"/tree_inference_depth_",max_depth,".stan")
  
  mod = cmdstan_model(file)
  
  fit = mod$variational(data = data, seed = seed,
                             init = init,
                             output_samples = ndraws, 
                             algorithm="fullrank")
  
  return(fit)
  
}


# Infer epimutation clocks in number of cell divisions and fitness of + cells with respect to - cells.
#
# A cmdstanr fit is returned.
#
# @param fit_tree cmdstanr fit containing tree inference
# @param path_to_model String specifying the path where we wannt to save the stan model for a given depth
# @param cmdstan_path String specifying the path to cmdstan folder
# @param threshold Threshold for tree pruning
# @param ndraws Number of draws from the posterior
# @param init List of initialization parameters
# @param seed Seed of the computation
# @param mu Mutation rate per division per bp per allele
# @param l length of the genome
#@ param ms Mean of lognormal prior for s
#@ param sigma Sigma parameter of lognormal prior for s
# @param k Number of trials for the beta prior on cluster centroids
# @return cdmstanr fit
# @examples
# fit_s(fit_tree,path_to_model = "models",cmdstan_path = "my_cmdstan/",threshold = 0.1,
# ndraws = 1000,init = NULL,seed = 45,
# mu = 1e-7,l = 2.7*10^9,ms = -0.5,sigma = 0.5,k = 100)
# @export

fit_s = function(fit_tree,path_to_model,cmdstan_path,threshold = 0.1,
                 ndraws = 1000,init = NULL,seed = 45,
                 mu = 1e-7,l = 2.7*10^9,ms = -0.5, sigma = 0.5, k = 100){
  
  dir.create(path_to_model)
  cmdstanr::set_cmdstan_path(cmdstan_path)
  
  tree = get_average_tree(fit_tree,threshold = threshold)
  
  model = fitness_inference_code(tree,likelihood = T)
  
  write_stan_file(
    model,
    dir = ".",
    basename = paste0(path_to_model,"/fitness_inference.stan"),
    force_overwrite = FALSE,
    hash_salt = ""
  )
  
  data = list(delta_t_n = 0,
              mu = mu,
              l = l,
              ms = ms,
              sigma = sigma,
              k = k)
  
  param = tree %>% reshape2::melt() %>% mutate(node = gsub(x = node,pattern = "-",replacement = "n")) %>% 
      mutate(node = gsub(x = node,pattern = "\\+",replacement = "p")) %>% filter(variable %in% c("m","vaf_minus","vaf_plus")) %>% 
      mutate(variable = paste0(variable,"_",node)) %>% dplyr::select(variable,value)
  
  extra_data = param$value
  names(extra_data) = param$variable
  
  data = append(data,extra_data)
  
  file = paste0(path_to_model,"/fitness_inference.stan")
  
  mod = cmdstan_model(file)
  
  fit = mod$variational(data = data, seed = seed,
                             init = init,
                             output_samples = ndraws, 
                             algorithm="fullrank")
  
  return(fit)
  
}


# x = fit_tree(spectrum,max_depth = 2,path_to_model = "models",cmdstan_path = "/opt/anaconda3/envs/stan/bin/cmdstan/",seed = 2000)
# y = fit_s(x,path_to_model = "models",cmdstan_path = "my_cmdstan/",threshold = 0.1,
# ndraws = 1000,init = NULL,seed = 45, mu = 1e-7,l = 2.7*10^9,ms = -0.5,sigma = 0.5,k = 100)
# tree = get_average_tree(x,threshold = NULL)
# clusters = get_clusters(spectrum %>% dplyr::select(-node),tree)




