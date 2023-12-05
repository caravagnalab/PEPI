# utilities inference

# Get the inferred tree from a fit.
#
# A tree is obtained from the inferred values of nodes ccfs and number of mutations.
#
# @param x A fit
# @param threshold Threshold on epimutation probability for tree pruning
# @return Inferred sample tree
# @examples
# get_average_tree(x,threshold = 0.1)
# @export


get_average_tree = function(x,threshold = NULL){
  
  params = x$summary() %>% as_tibble() %>% dplyr::rename(param = variable)
  M = round((x$draws("m_n") %>% as.vector())/(x$draws("nu_n") %>% as.vector())) %>% unique()
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
  
  return(tree)
  
}

# Get posterior draws from fit.
#
# Posterior draws are extracted from the fit.
#
# @param x Inference fit
# @param threshold Threshold for tree pruning
# @return A dataframe with posterior draws
# @examples
# get_posterior(x,threshold = 0.1)
# @export


get_posterior = function(x,model_type = "tree",threshold = 0.1){

if(model_type == "tree"){
  
  tree = get_average_tree(x,threshold = threshold)
  tree = tree %>% mutate(node = gsub(x = node,pattern = "-",replacement = "n")) %>% 
           mutate(node = gsub(x = node,pattern = "\\+",replacement = "p"))
  
  nodes = tree %>% pull(node)
  
 params =  tree %>% 
   mutate(w_plus  = ifelse(paste0(node,"-") %in% nodes,0.5,1), 
          w_minus = ifelse(paste0(node,"-") %in% nodes,0.5,1)) %>% reshape2::melt() %>% 
    filter(!(variable == "p_switch" & value == 0)) %>% 
    filter(!(variable == "nu" & value == 1)) %>% 
    filter(!(variable %in% c("w_minus","w_plus") & value == 1)) %>% 
    mutate(variable = paste0(variable,"_",node)) %>% 
     pull(variable)
 
   all_variables = x$draws() %>% as.data.frame()  %>%
                 reshape2::melt() %>% pull(variable) %>% unique()
   
   
  params = c(params[params %in% all_variables],"rn","rp")
  
  posterior = as.data.frame(x$draws())[,colnames(x$draws() %>% as.data.frame()) %in% params]

 full_tree =  get_average_tree(x,threshold = NULL)

if(nrow(full_tree) > nrow(tree)){
  
tree = tree %>% mutate(leave = ifelse(!paste0(node,"n") %in% nodes,T,F))
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
}
  
if(model_type == "fitness"){
  
  posterior =  x$draws() %>% as.data.frame() %>% as_tibble() %>% 
    dplyr::select(-c( "lp__","lp_approx__"))
  
}
  
  
return(posterior %>% as_tibble())
  
}

# Infer clusters from data and inferred tree.
#
# Membership probabilities and cluster assignents are computed for any mutation.
#
# @param data Dataset containing number of variants and depth for any mutation
# @param Inferred tree
# @return A list with cluster membership probabilities and mutation assignment
# @examples
# get_clusters(data,tree)
# @export

get_clusters = function(data,tree){
  
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
  
  return(list(ass_probs = probs, data = clusters))
  
}

# Get draws from prior used for the inference
#
# Prior draws are generated from inference object.
#
# @param x Fit obtained from the inference
# @param model_type: It can be "tree" or "fitness"
# @param threshold Threshold on switching probability for tree pruning
# @param ndraws Number of draws
# @return A dataframe with prior draws for any parameter
# @examples
# get_prior(x,model_type = "tree",ndraws = 1000)
# @export

# get_prior = function(x,model_type = "tree",threshold = 0.1,ndraws = 1000){
#   
#   data = fromJSON(file = x$data_file())
#   
#   if(model_type == "tree"){
#     max_depth = gsub(x = x$draws() %>% as.data.frame() %>% 
#                        dplyr::select(starts_with("vaf_minus")) %>% colnames(),
#                      pattern = "vaf_minus_",replacement = "") %>% nchar() %>% max() - 1
#     model_prior = tree_inference_code(max_depth,likelihood = F)
#   }
#    if(model_type == "fitness"){
#     tree = get_average_tree(x,threshold = threshold)
#     model_prior =  fitness_inference_code(tree,likelihood = F)
#   }
#   
#   prior <- rstan::stan(model_code = model_prior,data =  data, chains = 1, 
#                        warmup = 500, iter = ndraws,cores = 4)
# 
#   # model <- cmdstan_model(model_prior)
#   # prior <- model$sample(
#   #   data = data,
#   #   seed = 123,
#   #   chains = 1
#   # )
#   
#   prior = prior %>% as.data.frame() %>% as_tibble()
#   
#   return(prior)
#   
# }


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

plot_inference = function(post,params = NULL,prior = NULL){
  
  if(!is.null(prior)){ 
  prior = prior %>% as_tibble() %>% reshape2::melt() %>% 
          mutate(type = "prior")}
  if(!is.null(post)){
  post = post %>% as_tibble() %>% reshape2::melt() %>% 
               mutate(type = "post")}
  
  sampling = rbind(prior,post)  
  
  if(!is.null(params)){
    
    sampling = sampling %>% filter(variable %in% params) 
  }
  
  nr = round(length(sampling$variable %>% unique())/4 + 1)
  nc = min(length(sampling$variable %>% unique()) + 1,4) 
  
  ggplot(sampling) + geom_density(aes(x = value, alpha = type),fill = "steelblue") + 
     facet_wrap(~variable,scales = "free",nrow = nr, ncol = nc)  + 
    scale_alpha_manual(values = c("prior" = 0.4, "posterior" = 1)) + 
    CNAqc:::my_ggplot_theme() + theme(legend.position="none")
  
}









