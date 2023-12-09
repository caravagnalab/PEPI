

#' Fit a multivariate vaf spectrum with epigenetic tree model.
#'
#' A Pepi fit with tree inference is returned.
#'
#' @param x Pepi object containing VAF multivariate spectrum
#' @param path_to_model String specifying the path where we wannt to save the stan model for a given depth
#' @param cmdstan_path String specifying the path to cmdstan folder
#' @param max_depth Maximum number of levels
#' @param ndraws Number of draws from the posterior
#' @param init List of initialization parameters
#' @param seed Seed of the computation
#' @param mu Mutation rate per division per bp per allele
#' @param l length of the genome
#' @param rho_n Purity of - sample
#'  @param rho_p Purity of + sample
#' @param nu_t Mean of the beta prior on fraction of truncal mutations
#' @param qt Number of trials of the beta prior on fraction of truncal mutations
#' @param rate_n Mean of the beta prior on the epimutation rate from - to +
#' @param qn Number of trials for the beta prior on the epimutation rate from - to +
#' @param rate_p Mean of the beta prior on the epimutation rate from + to -
#' @param qp Number of trials for the beta prior on the epimutation rate from - to +
#' @param k Number of trials for the beta prior on cluster centroids
#' @param gamma Concentration of a Dirichlet distribution to split mutations at any node
#' @return Pepi object
#' @examples
#' fit_tree(x,path_to_model = "models",cmdstan_path = "my_cmdstan/",
#' max_depth = 2,ndraws = 1000,init = NULL,seed = 15,
#' mu = 1e-7,l = 2.7*10^9,rho_n = 1,rho_p = 1,nu_t = 0.1,
#' qt = 1e4,rate_n = 1e-3,qn = 1e4,rate_p = 1e-3,
#' qp = 1e4,k = 1e4,gamma = 150)
#' @export

fit_tree = function(x,path_to_model,cmdstan_path,
                    max_depth = 2,ndraws = 1000,init = NULL,seed = 15,
                    mu = 1e-7,l = 2.7*10^9,rho_n = 1,rho_p = 1,nu_t = 0.1,
                    qt = 1e4,rate_n = 1e-3,qn = 1e4,rate_p = 1e-3,
                    qp = 1e4,k = 1e4,gamma = 150){
  
  dir.create(path_to_model)
  cmdstanr::set_cmdstan_path(cmdstan_path)

  spectrum = x$VAF
  
  
  if(is.null(spectrum)){
    
    stop("no VAF spectrum") 
    
  }
  
  model = tree_inference_code(max_depth = max_depth,likelihood = T)

if(! paste0("tree_inference_depth_",max_depth,".stan") %in% list.files(path_to_model)){ 
  
  write_stan_file(
    model,
    dir = ".",
    basename = paste0(path_to_model,"/tree_inference_depth_",max_depth,".stan"),
    force_overwrite = FALSE,
    hash_salt = ""
  )
  
}
  
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
  
  pepi = list(inference = list(tree = fit),stan_data = list(tree = data), max_depth = max_depth)
  
  x$inference$tree = fit
  x$stan_data$tree = data
  x$max_depth = max_depth
  
  return(x)
  
}


#' Infer epimutation clocks in number of cell divisions and fitness of + cells with respect to - cells.
#'
#' A Pepi fit with fitness inference is returned.
#'
#' @param x Pepi object containing fitness and epimutation clocks inference
#' @param path_to_model String specifying the path where we want to save the stan model for a given depth
#' @param cmdstan_path String specifying the path to cmdstan folder
#' @param threshold Threshold for tree pruning
#' @param ndraws Number of draws from the posterior
#' @param init List of initialization parameters
#' @param seed Seed of the computation
#' @param mu Mutation rate per division per bp per allele
#' @param l length of the genome
#' @param ms Mean of lognormal prior for s
#'  @param sigma Sigma parameter of lognormal prior for s
#' @param k Number of trials for the beta prior on cluster centroids
#' @return Pepi object
#' @examples
#' fit_s(x,path_to_model = "models",cmdstan_path = "my_cmdstan/",threshold = 0.1,
#' ndraws = 1000,init = NULL,seed = 45,
#' mu = 1e-7,l = 2.7*10^9,ms = -0.5,sigma = 0.5,k = 100)
#' @export

fit_s = function(x,path_to_model,cmdstan_path,threshold = 0.1,
                 ndraws = 1000,init = NULL,seed = 45,
                 mu = 1e-7,l = 2.7*10^9,ms = -0.5, sigma = 0.5, k = 100){
  
  dir.create(path_to_model)
  cmdstanr::set_cmdstan_path(cmdstan_path)
  
  if(is.null(x$inference$tree)){
    
    stop("no tree inference") 
    
  }

  if(!is.null(x$inferred_tree)){  
  x = get_average_tree(x,threshold = threshold)}
  
  tree = x$inferred_tree
  
  model = fitness_inference_code(tree,likelihood = T)

if(! "/fitness_inference.stan" %in% list.files(path_to_model)){
  
  write_stan_file(
    model,
    dir = ".",
    basename = paste0(path_to_model,"/fitness_inference.stan"),
    force_overwrite = FALSE,
    hash_salt = ""
  )

}
  
  data = list(delta_t_n = 0,
              mu = mu,
              l = l,
              ms = ms,
              sigma = sigma,
              k = k)
  
  param = tree %>% reshape2::melt() %>% 
    mutate(node = gsub(x = node,pattern = "-",replacement = "n")) %>% 
      mutate(node = gsub(x = node,pattern = "\\+",replacement = "p")) %>% 
    filter(variable %in% c("m","vaf_minus","vaf_plus")) %>% 
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
  
  x$inference$fitness = fit
  x$stan_data$fitness = data
  
  return(x)
  
}




#' A Pepi fit with cell counts inference is returned.
#'
#' @param x Pepi object containing cell counts data
#' @param threshold Threshold for tree pruning
#' @param ndraws Number of draws from the posterior
#' @param init List of initialization parameters
#' @param seed Seed of the computation
#' @param alpha_ln shape parameter for gamma prior on - growth rate
#' @param beta_ln rate parameter for gamma prior on - growth rate
#' @param alpha_lp shape parameter for gamma prior on + growth rate
#' @param beta_lp rate parameter for gamma prior on + growth rate
#' @param alpha_rn shape parameter for gamma prior on - effective switch rate
#' @param beta_rn rate parameter for gamma prior on - effective switch rate
#' @param alpha_rp shape parameter for gamma prior on + effective switch rate
#' @param beta_rp rate parameter for gamma prior on + effective switch rate
#' @return Pepi fit
#' @examples
#' fit_counts(x,path_to_model,cmdstan_path,
#'                      ndraws = 1000,init = NULL,seed = 45, alpha_ln = 1.5,
#'                      beta_ln = 1, alpha_lp = 1.5, beta_lp = 1, alpha_rn = 1, beta_rn = 10,
#'                      alpha_rp = 1, beta_rp = 10)
#' @export

fit_counts = function(x,path_to_model,cmdstan_path,
                      ndraws = 1000,init = NULL,seed = 45, alpha_ln = 1.5,
                      beta_ln = 1, alpha_lp = 1.5, beta_lp = 1, alpha_rn = 1, beta_rn = 10,
                      alpha_rp = 1, beta_rp = 10
                      ){
  
  dir.create(path_to_model)
  cmdstanr::set_cmdstan_path(cmdstan_path)
  
  if(is.null(x$counts)){
    
    stop("no cell counts") 
    
  }
  
  model = counts_inference_code(likelihood = T)
  
  if(! "regressionODE.stan" %in% list.files(path_to_model)){
    
    write_stan_file(
      model,
      dir = ".",
      basename = paste0(path_to_model,"/regressionODE.stan"),
      force_overwrite = FALSE,
      hash_salt = ""
    )
    
  }
  
  t0 = counts$time %>% min()
  z0n = counts %>% filter(time == t0,epistate == "-") %>% pull(counts)
  z0p = counts %>% filter(time == t0,epistate == "+") %>% pull(counts)
    
  data  = list(
    n_times = counts %>% filter(time > t0) %>% pull(time) %>% unique() %>% length(),
    z0 = c(z0n,z0p,0,0,0),
    t0 = t0,
    zminus = counts %>% filter(time > t0, epistate == "-") %>% pull(counts), 
    zplus = counts %>% filter(time > t0, epistate == "+") %>% pull(counts),
    t = counts %>% filter(time > 0) %>% pull(time) %>% unique(),
    alpha_ln = alpha_ln,
    beta_ln = beta_ln,
    alpha_lp = alpha_lp,
    beta_lp = beta_lp,
    alpha_rn = alpha_rn,
    beta_rn = beta_rn,
    alpha_rp = alpha_rp,
    beta_rp = beta_rp 
)
  
 file = paste0(path_to_model,"/regressionODE.stan")
  
  mod = cmdstan_model(file)
  
  fit = mod$variational(data = data, seed = seed,
                        init = init,
                        output_samples = ndraws, 
                        algorithm="fullrank")
  
  x$inference$counts = fit
  x$stan_data$counts = data
  
  return(x)
  
}




