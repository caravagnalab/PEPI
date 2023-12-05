# generative models

# Generate a sample tree with input epimutation rates, fraction of truncal mutations and fitness advantage between + and -.
#
# A tree with number of mutations per nodes is generated.
#
# @param M total number of mutations
# @param gamma parameter of a beta distribution controlling how we split mutations at any node
# @param at number of successes of a beta distribution controlling amount of truncal mutations 
# @param bt number of fails of a beta distribution controlling amount of truncal mutations 
# @param s Fitness advantage of + cells with respect to - cells
# @param rate_minus Average number of epimutation from - to + per mutations
# @param rate_plus Average number of epimutation from + to - per mutations
# @param mu Mutation rate per division per bp per allele
# @param len length of the genome
# @param max_depth maximal level of the tree
# @return sample tree
# @examples
# generate_nodes(6000,100,10,80,0.2,1e-3,1e-3,1e-7,2.7*10^9,3)
# @export



generate_nodes = function(M,gamma,at,bt,s,rate_minus,rate_plus,mu,len,max_depth = 3){
  
  tree = data.frame(node = "-",level = 0, delta_m = M, phi = 1, p_switch = 1, nu = rbeta(1,at,bt), 
                    rate = rate_minus, s = 0) %>% 
    mutate(m = nu*M)
  
  for(l in 1:max_depth){
    
    epsilon = tree %>% filter(level == l-1) 
    node = epsilon %>% pull(node)
    muts = epsilon %>% pull(m)
    delta_muts = epsilon %>% pull(delta_m)
    
    
    new =  lapply(1:length(node), function(i){
      
      phi_n = rbeta(1,gamma,gamma)
       
      new_nodes = tibble(node = c(paste0(node[i],"-"),paste0(node[i],"+")), level = l,
                         phi = c(phi_n,1-phi_n),
                         rate = c(rate_minus,rate_plus), 
                         s= c(0,s)) %>% 
        mutate(delta_m = (delta_muts[i] - muts[i])*phi)
      
      if(l < max_depth){
        
        new_nodes =   new_nodes %>% rowwise() %>% 
          mutate(p_switch = 1 - exp(-rate*delta_m), 
                 nu = rbeta(n = 1, shape1 = (1-exp(-rate*delta_m))/rate, 
                            shape2 = delta_m - (1-exp(-rate*delta_m))/rate )) %>% 
          mutate(m = nu*delta_m)
        
      }else{
        
        new_nodes =  new_nodes %>% mutate(p_switch = 0, nu = 1, m = delta_m)   
        
      }
      
      new_nodes
      
  }) %>% bind_rows()
    
    tree = rbind(tree,new)
    
  }
  
  tree = tree %>% mutate(pi = m/M)
  
  tree = tree %>% rowwise() %>% mutate(delta_t = ifelse(level > 0,rgamma(1,shape = m, rate = 2*mu*len),0))
  
  return(tree)
  
}

# Prunes a sample tree.
#
# Nodes with low epimutation probability are removed and replaced with a leave.
#
# @param tree A sample tree 
# @param threshold Pruning threshold for epimutation probability
# @return sample tree
# @examples
# filter_nodes(tree,threshold = 0.01)
# @export

filter_nodes = function(tree,threshold = 0.01){
  
  max_depth = tree$level %>% max()
  
  nodes_list = c("-","--","-+")
  
  for(lev in 1:max_depth){
    
    nodes = tree %>% filter(level == max_depth - lev) %>% pull(node)
    
  for(nod in nodes){
      
      desc = c(paste0(nod,"-"),paste0(nod,"+"))
      
      p_switch = tree %>% filter(node == nod) %>% pull(p_switch)
      
    if(p_switch < threshold & ! sum(desc %in% nodes_list)){
        
        tree = tree %>% filter(! node %in% desc) %>% mutate(m = ifelse(node == nod,delta_m,m),
                                                            p_switch = ifelse(node == nod,0,p_switch),
                                                            nu =  ifelse(node == nod,1,nu))
    }else{
              nodes_list = c(nod,desc,nodes_list)
      }
    }
    
  }
  
  M = tree %>% pull(delta_m) %>% max()
  tree = tree %>% mutate(pi = m/M)
  return(tree)
  
}


# Calculate ccf for nodes of a sample tree.
#
# CCFs are computed according to number of mutations assigned to any node and fitness advantage.
#
# @param tree A sample tree 
# @param s Fitness advantage of + cells with respect to - cells
# @return sample tree
# @examples
# calculate_ccf(tree,s)
# @export

calculate_ccf = function(tree,s){
  
  if("ccf_plus" %in% colnames(tree) & "ccf_minus" %in% colnames(tree)){
    
    tree = tree %>% dplyr::select(-c("ccf_minus","ccf_plus","vaf_minus","vaf_plus"))  
  }
  
  max_depth = tree$level %>% max()
  all_nodes = tree %>% pull(node)
  leaves = tree %>% 
  mutate(leave = ifelse(!paste0(node,"+") %in% all_nodes & !paste0(node,"-") %in% all_nodes,T,F)) %>%
    filter(leave) %>% pull(node)
  coeff = list(1,1,1+s,1/(1+s))
  names(coeff) = c("--","++","-+","+-")
  
  ccf_leaves =  lapply(leaves,function(leave){
    
    t = tree %>% filter(node == leave) %>% pull(delta_t)
    last = substr(leave,start = nchar(leave),stop = nchar(leave))
    for (j in 1:(nchar(leave)-1)){
      ch = substr(leave,start = 1,stop = nchar(leave)-j)
      delta_t = tree %>% filter(node == ch) %>% pull(delta_t)
      t = t + delta_t*coeff[[paste0(substr(ch,start = nchar(ch),stop = nchar(ch)),last)]]
    }
    tibble(leave = leave,last = last, t = t, ccf = exp(-t))
  }) %>% bind_rows()
  
  tot = ccf_leaves %>% group_by(last) %>% summarize(n = sum(ccf))
  ccf_min = tot %>% filter(last == "-") %>% pull(n)
  ccf_pl = tot %>% filter(last == "+") %>% pull(n)
  ccf_leaves = ccf_leaves %>% mutate(ccf_minus = ifelse(last == "-",ccf/ccf_min,0)) %>% 
    mutate(ccf_plus = ifelse(last == "+",ccf/ccf_pl,0))
  
  ccf_nodes = full_join(tree,ccf_leaves %>% dplyr::rename(node = leave) %>% 
                          dplyr::select(node,ccf_minus,ccf_plus),by = "node")
  
  
  for(lev in 1:max_depth){
    
    nodes = ccf_nodes %>% filter(level == max_depth - lev,is.na(ccf_plus),is.na(ccf_minus)) %>% pull(node)
    
    for (nod in nodes){
      
      ccf_1 = c(ccf_nodes %>% filter(node == paste0(nod,"+")) %>% pull(ccf_plus),
                ccf_nodes %>% filter(node == paste0(nod,"+")) %>% pull(ccf_minus))
      ccf_2 = c(ccf_nodes %>% filter(node == paste0(nod,"-")) %>% pull(ccf_plus),
                ccf_nodes %>% filter(node == paste0(nod,"-")) %>% pull(ccf_minus))
      ccf = ccf_1 + ccf_2
      
      ccf_nodes  = ccf_nodes  %>% mutate(ccf_plus = ifelse(node == nod,ccf[1],ccf_plus),
                                         ccf_minus = ifelse(node == nod,ccf[2],ccf_minus))
    }
    
  }
  
  ccf_nodes = ccf_nodes %>% mutate(vaf_minus = ccf_minus/2,vaf_plus = ccf_plus/2)
  
  return(as_tibble(ccf_nodes))
  
}

# Associate number of tail mutations to the leaves of a sample tree.
#
# Number of tail mutations is assigned to any leave based on its ccf and mutation rate
#
# @param tree A sample tree 
# @param mu  Mutation rate
# @param len Length of the genome
# @param vaf_min Minimum detectable VAF of a mutation
# @return a table of mutations with number of variants and depth
# @examples
# gadd_tail_muts(tree,1e-07,2.7*10^9,0.05)
# @export

add_tail_muts = function(tree,mu,len = 2.7*10^9,vaf_min = 0.05){
  
  max_depth = tree %>% pull(level) %>% max()
  
  tree = tree %>% mutate(tail = 
                           ifelse(level == max_depth & vaf_minus > vaf_min, mu*len*(1/vaf_min - 1/vaf_minus),0)) %>% 
    mutate(tail = 
             ifelse(level == max_depth & vaf_plus > vaf_min, mu*len*(1/vaf_min - 1/vaf_plus),tail)) 
  
  return(tree)
  
}

# Simulate VAF multivariate spectrum associated to a sample tree.
#
# A bivariate binomial process is simulated any node, with number of trials and success probabilities given by number of mutations and ccfs. If tail mutations
# are included, thei VAF is sampled from a power law distribution.
#
# @param tree A sample tree 
# @param vaf_min Lower bound on VAF spectrum
# @param DP average coverage
# @return a table of mutations with number of variants and depth
# @examples
# generate_spectrum(tree,vaf_min = 0.05,DP = 150,tail = T)
# @export



generate_spectrum = function(tree,vaf_min,DP,tail = T){
  
  nodes = tree %>% pull(node)
  
  muts = lapply(nodes,function(nod){
    px = tree %>% filter(node == nod) %>% pull(vaf_minus)
    py = tree %>% filter(node == nod) %>% pull(vaf_plus)
    m = round(tree %>% filter(node == nod) %>% pull(m))
    DPx = rpois(m,lambda = DP)
    DPy = rpois(m,lambda = DP)
    if(px > 0){
      Nx = rbinom(m,size = DPx,prob = px)}else{ Nx = 0}
    if(py > 0){
      Ny = rbinom(m,size = DPy,prob = py)}else{ Ny = 0}
    
    tibble(node = nod,Nx = Nx,DPx = DPx,
           Ny = Ny,DPy = DPy)
    
  }) %>% bind_rows()
  
  if(tail){
    
  tail_nodes = tree %>% filter(tail > 0) %>% pull(node)
  
  tail_muts = lapply(tail_nodes,function(nod){
    px = tree %>% filter(node == nod) %>% pull(vaf_minus)
    py = tree %>% filter(node == nod) %>% pull(vaf_plus)
    m = round(tree %>% filter(node == nod) %>% pull(tail))
    DPx = rpois(m,lambda = DP)
    DPy = rpois(m,lambda = DP)
    if(px > 0){
      Nx = rbinom(m,size = DPx,prob = rtruncpareto(n = m,shape = 2,lower = vaf_min,upper = px))}else{ Nx = 0}
    if(py > 0){
      Ny = rbinom(m,size = DPy,prob = rtruncpareto(n = m,shape = 2,lower = vaf_min,upper = py))}else{ Ny = 0}
    
    tibble(node = paste0(nod,"tail"),Nx = Nx,Ny = Ny, DPx = DPx,DPy = DPy)
    
  }) %>% bind_rows()
  
    }else{
    
    tail_muts =  NULL
  }
  
  muts = rbind(muts,tail_muts)
  
  return(muts)
  
}


