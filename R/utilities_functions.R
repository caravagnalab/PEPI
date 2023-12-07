# utilities inference

# Get the inferred tree from a fit.
#
# A tree is obtained from the inferred values of nodes ccfs and number of mutations.
#
# @param x A pepi fit
# @param threshold Threshold on epimutation probability for tree pruning
# @return Pepi fit with inferred average sample tree
# @examples
# get_average_tree(x,threshold = 0.1)
# @export


get_average_tree = function(fit,threshold = NULL){
  
  x = fit$inferece$tree
  
  params = x$summary() %>% as_tibble() %>% dplyr::rename(param = variable)
  M = fit$stan_data$delta_m_n
  means = params %>% dplyr::select(param,mean)
  
  rn = means %>% filter(param == "rn") %>% pull(mean)
  rp = means %>% filter(param == "rp") %>% pull(mean)
  
  max_depth = gsub(x = params %>% filter(grepl(x = param,pattern = "vaf_minus")) %>% pull(param),
                pattern = "vaf_minus_",replacement = "") %>% 
                 nchar() %>% max() - 1
  
 tree = data.frame(node = "n",
                    level = 0, 
                    delta_m = M, 
                    phi = 1, 
                    p_switch = 1, 
                    nu = means %>% filter(param == "nu_n") %>% pull(mean),
                    rate = rn,
                    vaf_minus = means %>% filter(param == "vaf_minus_n") %>% pull(mean),
                    vaf_plus = means %>% filter(param == "vaf_plus_n") %>% pull(mean)) %>% mutate(m = nu*M)
  
  for(l in 1:max_depth){
    
    epsilon = tree %>% filter(level == l-1) 
    nod = epsilon %>% pull(node)
    muts = epsilon %>% pull(m)
    delta_muts = epsilon %>% pull(delta_m)
    
    
    new =  lapply(1:length(nod), function(i){
      
      phi =  means %>% filter(param == paste0("phi_",nod[i],"n")) %>% pull(mean)
      
      new_nodes = tibble(node = c(paste0(nod[i],"n"),paste0(nod[i],"p")), level = l,
                         rate = c(rn,rp), phi = c(phi,1-phi)) %>% rowwise() %>% 
        mutate(delta_m = (delta_muts[i] - muts[i])*phi, 
               vaf_minus = ifelse(paste0("vaf_minus_",node) %in% params$param,
                                  means %>% filter(param == paste0("vaf_minus_",node)) %>% pull(mean),0),
               vaf_plus = ifelse(paste0("vaf_plus_",node) %in% params$param,
                                 means %>% filter(param == paste0("vaf_plus_",node)) %>% pull(mean),0))
      
      if(l < max_depth){
        
        new_nodes = new_nodes %>% rowwise() %>% 
          mutate(p_switch = 1 - exp(-rate*delta_m), 
                 nu = means %>% filter(param == paste0("nu_",node)) %>% pull(mean)
          ) %>% 
          mutate(m = nu*delta_m)
        
      }else{
        
        new_nodes =  new_nodes %>% 
          mutate(p_switch = 0, nu = 1, m = delta_m)
        
      }
      new_nodes 
      
    }) %>% bind_rows()
    
    tree = rbind(tree,new)
    
  }
  
  tree$node = gsub(x = tree$node,pattern = "n",replacement = "-")
  tree$node = gsub(x = tree$node,pattern = "p",replacement = "+")
  tree = tree %>% mutate(pi = m/M) %>% as_tibble()
  
  # tree pruning
  
if(!is.null(threshold)){
  
  tree = filter_nodes(tree,threshold) 
  nodes = tree %>% pull(node)
  tree = tree %>% mutate(leave = ifelse(!paste0(node,"-") %in% nodes,T,F))
  leaves = tree %>% filter(leave) %>% pull(node)
  
  for (leav in leaves){
    
    v_m = ifelse(substr(x = leav,start = nchar(leav),stop = nchar(leav)) == "-",
                 tree %>% filter(node == substr(x = leav,start = 1,stop = nchar(leav)-1)) %>% 
                   pull(vaf_minus),1e-6)
    v_p = ifelse(substr(x = leav,start = nchar(leav),stop = nchar(leav)) == "+",
                 tree %>% filter(node == substr(x = leav,start = 1,stop = nchar(leav)-1)) %>% 
                   pull(vaf_plus),1e-6)
    
    delt_m = (tree %>% filter(node == substr(x = leav,start = 1,stop = nchar(leav)-1)) %>% 
                pull(delta_m) - tree %>% filter(node == substr(x = leav,start = 1,stop = nchar(leav)-1)) %>% 
                 pull(m))*(tree %>% filter(node == leav) %>% pull(phi))
      
    tree = tree %>% rowwise() %>% 
        mutate(vaf_minus = ifelse(node == leav,v_m,vaf_minus),
               vaf_plus = ifelse(node == leav,v_p,vaf_plus),
               delta_m = ifelse(node == leav,delt_m,delta_m),
               m = ifelse(node == leav,delt_m,m)) %>% ungroup()

}    
  }
  
  fit$inferred_tree = tree
  
  fit$epimutation_threshold = threshold
  
  return(fit)
  
}

# Get posterior draws from fit.
#
# Posterior draws are extracted from the fit.
# @param threshold Threshold on epimutation probability for tree pruning
# @param x Pepi fit
# @return A pepi object with posterior draws
# @examples
# get_posterior(fit)
# @export


get_average_counts = function(x){
  
  if(is.null(x$inference$counts){
    
    stop("no inference counts")
    
  }
  
  x$inference$counts$draws() %>% as.data.frame() %>% dplyr::select(starts_with("pred_minus")) %>% 
    reshape2::melt() %>% group_by(variable) %>% summarize(counts = mean(value))  %>% mutate(time = x$stan_data$counts$t)
  
}

get_posterior = function(fit,threshold = NULL){

if("tree" %in% names(x$inference)){
  
  if("inferred_tree" %in% names(fit)){
    
    fit = get_average_tree(fit,threshold)
    
  }
  
  x = fit$inference$tree
  
  tree = fit$inferred_tree
  
  tree = tree %>% mutate(node = gsub(x = node,pattern = "-",replacement = "n")) %>% 
           mutate(node = gsub(x = node,pattern = "\\+",replacement = "p"))
  
  tree = tree %>% mutate(leave = ifelse(!paste0(node,"n") %in% nodes,T,F))
  
  nodes = tree %>% pull(node)
  
 params =  tree %>% 
   mutate(w_plus  = 0.5, w_minus = 0.5) %>% reshape2::melt() %>% 
    filter(!(variable == "p_switch" & leave)) %>% 
    filter(!(variable == "nu" & leave)) %>% 
    filter(!(variable %in% c("w_minus","w_plus") & leave)) %>% 
    filter(!(variable %in% c("w_minus","w_plus") & 
    substr(x = node,start = nchar(node),stop = nchar(node)) == "p"),
     !(node == "n" & variable %in% c("phi","w_minus","w_plus","p_switch")),
     ! variable %in% c("level","s","delta_t","rate","tail","ccf_minus","ccf_plus")) %>% 
    mutate(variable = paste0(variable,"_",node)) %>% pull(variable)
 
  params = c(params,"rn","rp")
    
  posterior = as.data.frame(x$draws())[,colnames(x$draws() %>% as.data.frame()) %in% params]

  max_depth = fit$max_depth 

if( nrow(tree) <  2**(max_depth+1) - 1){
  
leaves = tree %>% filter(leave) %>% pull(node)
  
for (leav in leaves){

posterior[,colnames(posterior) == paste0("delta_m_",leav)] = 
    (posterior[,colnames(posterior) == paste0("delta_m_",substr(x = leav,start = 1,stop = nchar(leav)-1))] - 
       posterior[,colnames(posterior) == paste0("m_",substr(x = leav,start = 1,stop = nchar(leav)-1))])*(
         posterior[,colnames(posterior) == paste0("phi_",substr(x = leav,start = 1,stop = nchar(leav)))])
     
   posterior[,colnames(posterior) == paste0("m_",leav)] = 
      posterior[,colnames(posterior) == paste0("delta_m_",leav)]
    
    if(substr(x = leav,start = 1,stop = nchar(leav)) == "n"){
    posterior[,colnames(posterior) == paste0("vaf_minus_",leav)] = 
    posterior[,colnames(posterior) == paste0("vaf_minus_",substr(x = leav,start = 1,stop = nchar(leav)-1))]
      }else{
    posterior[,colnames(posterior) == paste0("vaf_minus_",leav)] =  rbeta(nrow(posterior),1,1e6)
      }
    
    if(substr(x = leav,start = 1,stop = nchar(leav)) == "p"){
      posterior[,colnames(posterior) == paste0("vaf_plus_",leav)] = 
        posterior[,colnames(posterior) == paste0("vaf_plus_",substr(x = leav,start = 1,stop = nchar(leav)-1))]
    }else{
      posterior[,colnames(posterior) == paste0("vaf_plus_",leav)] =  rbeta(nrow(posterior),1,1e6)
    }

    }    
}
  
   fit$posterior$tree = posterior
}
  
if("fitness" %in% names(x$inference)){
  
  x = fit$inference$fitness
  posterior =  x$draws() %>% as.data.frame() %>% as_tibble() %>% 
    dplyr::select(-c( "lp__","lp_approx__"))
  
  fit$posterior$fitness = posterior
  
}
  
if("counts" %in% names(x$inference)){
    
  x = fit$inference$counts
    posterior =  x$draws() %>% as.data.frame() %>% as_tibble() %>% 
      dplyr::select(-c( "lp__","lp_approx__"))
    
  fit$posterior$counts = posterior
    
  }
  
 return(fit)
  
}

# Infer clusters from data and inferred tree.
#
# Membership probabilities and cluster assignents are computed for any mutation.
#
# @param fit Pepi object
# @return A pepi object containing membership probabilities and cluster assignents
# @examples
# get_clusters(fit)
# @export

get_clusters = function(fit){
  
  data = fit$VAF
  tree = fit$inferred_tree
  
  if(is.null(tree)){
    
    stop("no inferred tree")
    
  }
  
  data = data %>% mutate(VAFx = Nx/DPx + 1e-3,VAFy = Ny/DPy + 1e-3) %>% 
    mutate(data_id = 1:nrow(data))
  
  if("node" %in% colnames(data)){
    
    data = data %>% dplyr::select(-node)
    
  }
  
  tree = tree %>% mutate(vaf_plus = vaf_plus + 1e-3,vaf_minus = vaf_minus + 1e-3)
  
  ass = lapply(1:nrow(tree), function(i){
    
    llk = data %>% mutate(llk = log(tree[i,]$pi) + 
                            dbeta(VAFx,shape1 = tree[i,]$vaf_minus*DPx,
                                  shape2 = (1-tree[i,]$vaf_minus)*DPx,
                                  log = T) + 
                            dbeta(VAFy,shape1 = tree[i,]$vaf_plus*DPy,
                                  shape2 = (1-tree[i,]$vaf_plus)*DPy,
                                  log = T)) %>% 
      mutate(node = tree$node[i])
    
    
    llk %>%  dplyr::select(llk,node,data_id)
    
  }) %>% bind_rows()
  
  probs = mclapply(1:nrow(data),function(i){
    
    x = ass %>% filter(data_id == i)
    
    tot = max(x$llk) + log(sum(exp(x$llk-max(x$llk))))
    
    x$llk = exp(x$llk - tot)
    
    x  
    
  }) %>% bind_rows()
  
  clusters = mclapply(1:nrow(data),function(i){
    
    x = probs %>% filter(data_id == i)
    
    cl =  x[x$llk == max(x$llk),]$node
    
    tibble(data_id = i, node= cl)
    
  }) %>% bind_rows()
  
  clusters = full_join(data,clusters,by = "data_id")
  
  fit$cluster_probs = probs
  
  fit$VAF = clusters
    
  return(fit)
  
}

# to fix

# Get draws from prior used for the inference
#
# Prior draws are generated from Pepi object.
#
# @param Pepi object
# @param model_types list of models we want the prior
# @param ndraws Number of draws
# @return A Pepi object with prior draws
# @examples
# get_prior(x,model_type = c("tree","counts"),ndraws = 1000)
# @export

get_prior = function(x,model_types = c("tree","counts"),ndraws = 1000){

 if("counts" %in% model_types & ! is.null(x$stan_data$counts)){
    
 model_prior =  counts_inference_code(likelihood = F)
  sampling <- rstan::stan(model_code = model_prior,data =  x$stan_data$counts, chains = 1,
                       warmup = 500, iter = ndraws,cores = 4)
  x$prior$counts = sampling %>% as.data.frame() %>% as_tibble() %>% 
    dplyr::select(-c( "lp__"))
 }

  if("tree" %in% model_types & ! is.null(x$stan_data$tree) ){
    
    max_depth = x$max_depth
    model_prior = tree_inference_code(max_depth,likelihood = F)
    sampling <- rstan::stan(model_code = model_prior,data =  x$stan_data$tree, chains = 1,
                            warmup = 500, iter = ndraws,cores = 4)
    x$prior$tree = sampling %>% as.data.frame() %>% as_tibble() %>% 
      dplyr::select(-c( "lp__"))
  }
   if("fitness" %in% model_types & ! is.null(x$stan_data$fitness)){
     
    tree = x$inferred_tree
    model_prior =  fitness_inference_code(tree,likelihood = F)
    sampling <- rstan::stan(model_code = model_prior,data =  x$stan_data$fitness, chains = 1,
                            warmup = 500, iter = ndraws,cores = 4)
    x$prior$fitness = sampling %>% as.data.frame() %>% as_tibble() %>% 
      dplyr::select(-c( "lp__"))
  }

    return(x)

}


# Plot posterior and prior distributions.
#
# Posterior and prior draws histograms are plotted for any required parameter.
#
# @param post Posterior draws
# @param params A list of parameters 
# @return A plot with posterior and prior distributions
# @examples
# plot_inference(x,params = c("rn","rp"))
# @export

plot_inference = function(x,params = NULL){
  
  prior = rbind(as_tibble(x$prior$tree) %>% mutate(class = "tree"),
                as_tibble(x$prior$fitness) %>% mutate(class = "fitness"),
                as_tibble(x$prior$counts) %>% mutate(class = "counts")) %>% 
              mutate(type = "prior")  %>% reshape2::melt()
  
  post = rbind(as_tibble(x$post$tree) %>% mutate(class = "tree"),
                as_tibble(x$post$fitness) %>% mutate(class = "fitness"),
                as_tibble(x$post$counts) %>% mutate(class = "counts")) %>% 
              mutate(type = "post")  %>% reshape2::melt()
  
  sampling = rbind(prior,post)
  
  if(!is.null(params)){
    
    sampling = sampling %>% filter(variable %in% params) 
  }
  
  if (is.null(sampling)) {
    stop("required parameters are not present")
    
  }
  
  nr = round(length(sampling$variable %>% unique())/4 + 1)
  nc = min(length(sampling$variable %>% unique()) + 1,4) 
  cls = c("indianred","steelblue","darkgreen")
  names(cls) = c("tree","fitness","counts")
    
  ggplot(sampling) + geom_density(aes(x = value, alpha = type,fill = class)) + 
     facet_wrap(~variable,scales = "free",nrow = nr, ncol = nc)  + 
    scale_alpha_manual(values = c("prior" = 0.4, "posterior" = 1)) + 
    scale_fill_manual(values = cls) + 
    CNAqc:::my_ggplot_theme() + theme(legend.position="none")
  
}









