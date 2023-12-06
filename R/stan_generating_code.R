
# Generate stan model for tree inference.
#
# A stan file with the tree inference model is generated.
#
# @param max_depth Depth of the tree 
# @param likelihood Boolean value specifying to include the likelihood in the model, or just sampling from the priors.
# @return A string containing a stan model
# @examples
# tree_inference_code(max_depth = 2,likelihood = T)
# @export



tree_inference_code = function(max_depth,likelihood = T){
  
  tree = tibble(node = "n",level = 0, delta_m = "delta_m_n", phi = "", p_switch = "", nu = "nu_n") %>% 
    mutate(m = "m_n = nu_n*delta_m_n", sampling_nu = "nu_n ~ beta_proportion(nu_t,qt)", sampling_phi = "")
  
  for(l in 1:max_depth){
    
    nod = tree %>% filter(level == l-1) %>% pull(node)
    
    new =  lapply(1:length(nod), function(i){
      
      new_nodes = tibble(node = c(paste0(nod[i],"n"),paste0(nod[i],"p")), level = l,
                         rate = c("rn","rp"), 
                         phi = c(paste0("phi_",nod[i],"n"),paste0("phi_",nod[i],"p")),
                         sampling_phi = 
                           c(paste0("phi_",nod[i],"n"," ~ beta(gamma,gamma)"),
                             paste0("phi_",nod[i],"p = 1 - ","phi_",nod[i],"n"))) %>% 
        mutate(delta_m = paste0("delta_m_",node,"=(delta_m_",nod[i],"- m_",nod[i],")*",phi), 
               nu = paste0("nu_",node))
      
      if(l < max_depth){
        
        new_nodes =   new_nodes %>% rowwise() %>% 
          mutate(p_switch = paste0("p_switch_",node," = ","1-exp(-",rate,"*delta_m_",node,")"),
                 sampling_nu = paste0("nu_",
                                      node," ~ beta(p_switch_",node,"/",rate,",delta_m_",node,"- p_switch_",node,"/",rate,")")) %>% 
          mutate(m = paste0("m_",node,"= ","nu_",node,"*","delta_m_",node))
        
      }else{
        
        new_nodes =  new_nodes %>% 
          mutate(p_switch = "", sampling_nu = "", m = paste0("m_",node,"= delta_m_",node))   
        
      }
      
      new_nodes = new_nodes %>% dplyr::select(-rate)
      
    }) %>% bind_rows()
    
    tree = rbind(tree,new)
    
  }
 
# vaf statement 
   
 
  tree = tree %>% mutate(vaf_minus = paste0("real <lower = 0,upper = 1> vaf_minus_",node,";"),
                         vaf_plus = paste0("real <lower = 0,upper = 1> vaf_plus_",node,";"))
  
  tree = tree %>% mutate(
         vaf_minus = ifelse(level < max_depth & level > 0 &
                           substr(node,start = nchar(node),stop = nchar(node)) == "n",
                           paste0("real <lower = 0,upper = ","vaf_minus_",
                 substr(node,start = 1,stop = nchar(node)-1),"> vaf_minus_",
                node," = w_minus_",substr(node,start = 1,stop = nchar(node)-1),"n*vaf_minus_",substr(node,start = 1,stop = nchar(node)-1),";"),vaf_minus),
          vaf_minus = ifelse(level < max_depth & level > 0 &
                               substr(node,start = nchar(node),stop = nchar(node)) == "p",
                             paste0("real <lower = 0,upper = ","vaf_minus_",
                                    substr(node,start = 1,stop = nchar(node)-1),"> vaf_minus_",node," = (1-w_minus_",
                                    substr(node,start = 1,stop = nchar(node)-1),"n)*vaf_minus_",
                                    substr(node,start = 1,stop = nchar(node)-1),";"),vaf_minus),
         vaf_plus = ifelse(level < max_depth & level > 0 &
                              substr(node,start = nchar(node),stop = nchar(node)) == "n",
                            paste0("real <lower = 0,upper = ","vaf_plus_",
                                   substr(node,start = 1,stop = nchar(node)-1),"> vaf_plus_",
                                   node," = w_plus_",substr(node,start = 1,stop = nchar(node)-1),"n*vaf_plus_",substr(node,start = 1,stop = nchar(node)-1),";"),vaf_plus),
         vaf_plus = ifelse(level < max_depth & level > 0 &
                             substr(node,start = nchar(node),stop = nchar(node)) == "p",
                           paste0("real <lower = 0,upper = ","vaf_plus_",
                                  substr(node,start = 1,stop = nchar(node)-1),"> vaf_plus_",node," = (1-w_plus_",
                                  substr(node,start = 1,stop = nchar(node)-1),"n)*vaf_plus_",
                                  substr(node,start = 1,stop = nchar(node)-1),";"),vaf_plus)
        
          )
  

  tree = tree %>% mutate(vaf_minus_sampling = 
                        ifelse(node == "n",paste0("vaf_minus_",node,"~ beta_proportion(0.5*rho_n,k);"),""),
                         vaf_plus_sampling = 
                        ifelse(node == "n",paste0("vaf_plus_",node,"~ beta_proportion(0.5*rho_p,k);"),""),
                        w_minus = "",
                        w_plus = "",
                        w_minus_sampling = "",
                        w_plus_sampling = "")
  
for(lev in 0:(max_depth - 1)){
  
  nodes = tree %>% filter(level == lev) %>% pull(node)
  
for (j in 1:length(nodes)){

    if(lev + 1 < max_depth){   
      
      tree =  tree %>% mutate(w_minus = ifelse(node == paste0(nodes[j],"n"),
                                               paste0("real <lower = 0,upper = 1 > ","w_minus_",nodes[j],"n;"),
                                               w_minus),
                              w_plus = ifelse(node == paste0(nodes[j],"n"),
                                               paste0("real <lower = 0,upper = 1 > ","w_plus_",nodes[j],"n;"),
                                               w_plus),
                              w_minus_sampling = ifelse(node == paste0(nodes[j],"n"),
                                                          paste0("w_minus_",paste0(nodes[j],"n"),"~ uniform(0,1);"),
                                                          w_minus_sampling),
                              w_plus_sampling = ifelse(node == paste0(nodes[j],"n"),
                                                         paste0("w_plus_",paste0(nodes[j],"n"),"~ uniform(0,1);"),
                                                         w_plus_sampling))
    }else{
    
  tree =  tree %>% mutate(vaf_minus_sampling = ifelse(node == paste0(nodes[j],"n"),
              paste0("vaf_minus_",paste0(nodes[j],"n"),
              "~ beta_proportion(vaf_minus_",nodes[j]," + 1e-3,k);"),vaf_minus_sampling),
                          vaf_plus_sampling = ifelse(node == paste0(nodes[j],"n"),
        paste0("vaf_plus_",paste0(nodes[j],"n"),"~ beta(1,1e6);"),vaf_plus_sampling))
  
  tree =  tree %>% mutate(vaf_plus_sampling = ifelse(node == paste0(nodes[j],"p"),
                  paste0("vaf_plus_",paste0(nodes[j],"p"),
                         "~ beta_proportion(vaf_plus_",nodes[j]," + 1e-3,k);"),
                  vaf_plus_sampling),
                    vaf_minus_sampling = ifelse(node == paste0(nodes[j],"p"),
                 paste0("vaf_minus_",paste0(nodes[j],"p"),"~ beta(1,1e6);"),vaf_minus_sampling))
  
       }
  
    }
}
  
  # data
  
  stick_breaking = 
    "
data{

real <lower=0> delta_m_n;
real <lower=0,upper = 1> nu_t;
real<lower=0> qt;
real <lower=0, upper = 1> rate_n;
real<lower=0> qn;
real<lower=0, upper = 1> rate_p;
real<lower=0> qp;
int n;
array[n] int Nn;
array[n] int Np;
array[n] int DPn;
array[n] int DPp;
real <lower=0> gamma;
real <lower=0> k;
real <lower=0> mu;
real <lower=0> l;
real <lower=0,upper = 1> rho_n;
real <lower=0,upper = 1> rho_p;

" 
  
  # params
  
  stick_breaking = paste0(stick_breaking, 
                          " }
                        
                        
parameters{
real <lower=0,upper = 1> rn;
real <lower=0,upper = 1 > rp;
  "
  )
  
  for (i in 1:nrow(tree)){

    if(tree$vaf_minus_sampling[i] != ""){    
    stick_breaking = paste(stick_breaking,tree$vaf_minus[i],"\n")}
    if(tree$vaf_plus_sampling[i] != ""){
    stick_breaking = paste(stick_breaking,tree$vaf_plus[i],"\n")}
    
    if(tree$w_minus_sampling[i] != ""){    
      stick_breaking = paste(stick_breaking,tree$w_minus[i],"\n")}
    if(tree$w_plus_sampling[i] != ""){
      stick_breaking = paste(stick_breaking,tree$w_plus[i],"\n")}
    
    if((i %% 2) == 0){
      stick_breaking = paste(stick_breaking,"real <lower=0, upper=1> ",tree$phi[i],";\n")}
    
    if(tree$level[i] < max_depth){
      stick_breaking = paste(stick_breaking,"real <lower=0, upper=1> ",tree$nu[i],";\n")}
    
  }
  
  # transformed params
  
  stick_breaking = paste0(stick_breaking," } 

transformed parameters{ \n")
  
  for (i in 1:nrow(tree)){
    
    if(tree$vaf_minus_sampling[i] == ""){    
      stick_breaking = paste(stick_breaking,tree$vaf_minus[i],"\n")}
    if(tree$vaf_plus_sampling[i] == ""){
      stick_breaking = paste(stick_breaking,tree$vaf_plus[i],"\n")}
    
    if((i %% 2) == 1 & i>1){
      stick_breaking = paste(stick_breaking,"real <lower=0, upper= 1>",tree$sampling_phi[i],";\n")}
    
    if( i > 1){
      
      stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$delta_m[i],";\n")}
    
    stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$m[i],";\n")
    
    if(i > 1 & tree$level[i] < max_depth){
      stick_breaking = paste0(stick_breaking,"real <lower=0, upper= 1>",tree$p_switch[i],";\n")}
    
  }
  
  stick_breaking = paste0(stick_breaking, "

} 

model{
rn ~ beta_proportion(rate_n,qn);
rp ~ beta_proportion(rate_p,qp); \n")
  
  for (i in 1:nrow(tree)){
    
    if(tree$w_minus_sampling[i] != ""){    
      stick_breaking = paste(stick_breaking,tree$w_minus_sampling[i] ,"\n")}
    if(tree$w_plus_sampling[i] != ""){
      stick_breaking = paste(stick_breaking,tree$w_plus_sampling[i] ,"\n")}

    if(tree$vaf_minus_sampling[i] != ""){    
    stick_breaking = paste(stick_breaking,tree$vaf_minus_sampling[i] ,"\n")}
    if(tree$vaf_plus_sampling[i] != ""){
    stick_breaking = paste(stick_breaking,tree$vaf_plus_sampling[i] ,"\n")}
    
    if((i %% 2) == 0){
      stick_breaking = paste(stick_breaking,tree$sampling_phi[i],";\n")
    }
    
    if(tree$level[i] < max_depth){
      stick_breaking = paste0(stick_breaking,tree$sampling_nu[i],";\n")
    }
    
  }

if(likelihood){
  
  stick_breaking = paste0(stick_breaking,"

for (i in 1:n){

target += log_sum_exp([log(m_n/delta_m_n) + binomial_lpmf(Nn[i] | DPn[i], vaf_minus_n) + 
                                            binomial_lpmf(Np[i] | DPp[i], vaf_plus_n)\n ")
  
  for (i in 2:nrow(tree)){
    
    stick_breaking = paste0(stick_breaking,"        , log(m_",tree$node[i],"/delta_m_n) + \n  ",
                            "        binomial_lpmf(Nn[i] | DPn[i], vaf_minus_",tree$node[i],") + \n",
                            "        binomial_lpmf(Np[i] | DPp[i], vaf_plus_",tree$node[i],")\n")
    
    
  }
  
  
  # generated quantities
  
  stick_breaking = paste0(stick_breaking,
                          " 
]);
   }")
  
}  
  
stick_breaking = paste0(stick_breaking,"

}

generated quantities{
                        
     real<lower=0> effective_switch_rate_n = rn*2*mu*l;
     real<lower=0> effective_switch_rate_p = rp*2*mu*l;                   
                        
    }")
  
  return(stick_breaking)
  
}


# Generate stan model for fitness inference.
#
# A stan file with the fitness inference model is generated.
#
# @param tree Sample tree with epigenetic state associated to any node 
# @param likelihood Boolean value specifying to include the likelihood in the model, or just sampling from the priors.
# @return A string containing a stan model
# @examples
# fitness_inference_code(tree,likelihood = T)
# @export

fitness_inference_code = function(tree,likelihood = T){
  
  tree = tree %>% mutate(node = gsub(x = node,pattern = "-","n")) %>% 
    mutate(node = gsub(x = node,pattern = "\\+","p"))
  
  
  tree = tree %>% rowwise() %>% mutate(delta_t = paste0("delta_t_",node), 
                                       sampling_delta_t = paste0(delta_t, " ~ gamma(m_",node,",2*mu*l)"))
  
  max_deepth = tree$level %>% max()
  all_nodes = tree %>% pull(node)
  leaves = tree %>% mutate(leave = ifelse(!paste0(node,"p") %in% all_nodes & !paste0(node,"n") %in% all_nodes,T,F)) %>%
    filter(leave) %>% pull(node)
  coeff = list("1","1","(1+s)","1/(1+s)")
  names(coeff) = c("nn","pp","np","pn")
  
  ccf_leaves =  lapply(leaves,function(leave){
    
    t = tree %>% filter(node == leave) %>% pull(delta_t)
    last = substr(leave,start = nchar(leave),stop = nchar(leave))
    for (j in 1:(nchar(leave)-1)){
      ch = substr(leave,start = 1,stop = nchar(leave)-j)
      delta_t = tree %>% filter(node == ch) %>% pull(delta_t)
      t = paste0(t,"+",delta_t,"*",coeff[[paste0(substr(ch,start = nchar(ch),stop = nchar(ch)),last)]])
    }
    tibble(leave = leave,last = last, t = t, ccf = paste0("exp(-(",t,"))"))
  }) %>% bind_rows()
  
  tot = ccf_leaves %>% group_by(last) %>% summarize(n = paste0(ccf,collapse = "+"))
  ccf_minus = tot %>% filter(last == "n") %>% pull(n)
  ccf_plus = tot %>% filter(last == "p") %>% pull(n)
  ccf_leaves = ccf_leaves %>% mutate(ccf_minus = ifelse(last == "n",
                                                        paste0(ccf,"/(",ccf_minus,")"),paste0("2*vaf_minus_",leave))) %>% 
    mutate(ccf_plus = ifelse(last == "p",paste0(ccf,"/(",ccf_plus,")"),paste0("2*vaf_plus_",leave)))
  
  if("ccf_minus" %in% colnames(tree) | "ccf_plus" %in% colnames(tree)){
     tree = tree %>% dplyr::select(-starts_with("ccf"))
  }
    
  tree = full_join(tree,ccf_leaves %>% dplyr::rename(node = leave) %>% 
                     dplyr::select(node,ccf_minus,ccf_plus),by = "node")
  
  for(lev in 1:max_deepth){
    
    nodes = tree %>% filter(level == max_deepth - lev,
                            is.na(ccf_plus),
                            is.na(ccf_minus)) %>% pull(node)
    
    for (nod in nodes){
      
      # ccf_1 = c(tree %>% filter(node == paste0(nod,"p")) %>% pull(ccf_plus),
      #           tree %>% filter(node == paste0(nod,"p")) %>% pull(ccf_minus))
      # ccf_2 = c(tree %>% filter(node == paste0(nod,"n")) %>% pull(ccf_plus),
      #           tree %>% filter(node == paste0(nod,"n")) %>% pull(ccf_minus))
      # 
      # ccf = paste0(ccf_1,"+",ccf_2)
      
      ccf_1 = c(paste0("ccf_plus_",paste0(nod,"p")),paste0("ccf_minus_",paste0(nod,"p")))
      ccf_2 = c(paste0("ccf_plus_",paste0(nod,"n")),paste0("ccf_minus_",paste0(nod,"n")))
      ccf = paste0(ccf_1,"+",ccf_2)
      
      tree  = tree  %>% mutate(ccf_plus = ifelse(node == nod,ccf[1],ccf_plus),
                               ccf_minus = ifelse(node == nod,ccf[2],ccf_minus))
    }
    
  }
  
  tree = tree %>% mutate(ccf_minus = paste0("ccf_minus_",node,"=",ccf_minus),
                         ccf_plus = paste0("ccf_plus_",node,"=",ccf_plus))
  
  
  # data
  
stick_breaking = 
    "
data{

real <lower=0> delta_t_n;
real l;
real mu;
real ms;
real <lower=0> sigma;
real <lower=0> k;
" 
  
  for (i in 1:nrow(tree)){
    
    stick_breaking = paste0(stick_breaking,"real <lower=0>",
                            "m_",tree$node[i],";\n")
    
    stick_breaking = paste0(stick_breaking,"real <lower=0,upper = 1>",
                              "vaf_minus_",tree$node[i],";\n")
    
    
    stick_breaking = paste0(stick_breaking,"real <lower=0,upper = 1>",
                              "vaf_plus_",tree$node[i],";\n")
    
}
  
  # params
  
  stick_breaking = paste0(stick_breaking, 
                          " }
                        
                        
parameters{

real <lower=0> s;
  "
  )
  
  for (i in 2:nrow(tree)){
    
    stick_breaking = paste(stick_breaking,"real <lower=0> ",tree$delta_t[i],";\n")
    
  }
  
  # transformed params
  
  stick_breaking = paste0(stick_breaking," } 

transformed parameters{ 
                        ")
  
  for (i in 1:nrow(tree)){
    
    k = nrow(tree) - i + 1
    
    stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$ccf_minus[k],";\n")
      
    
    stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$ccf_plus[k],";\n")
 
  }
  
  stick_breaking = paste0(stick_breaking, "

} 

model{
       s ~  lognormal(ms,sigma); \n")
  
  for (i in 2:nrow(tree)){
    
    
    stick_breaking = paste0(stick_breaking,tree$sampling_delta_t[i],";\n")
    
  }
  

if(likelihood){
  
  for (i in 1:nrow(tree)){
 
  if(!grepl(x = tree$ccf_minus[i],pattern = "vaf")){  
    stick_breaking = paste0(stick_breaking," vaf_minus_",tree$node[i],
                            " ~ beta_proportion(0.5*ccf_minus_",tree$node[i],",k);\n")}
  
    if(!grepl(x = tree$ccf_plus[i],pattern = "vaf")){ 
    stick_breaking = paste0(stick_breaking," vaf_plus_",tree$node[i],
                            " ~ beta_proportion(0.5*ccf_plus_",tree$node[i],",k);\n")}
    
  }

}
  
  stick_breaking = paste0(stick_breaking,"
                        
 } ")
  
  return(stick_breaking)
  
}


# Generate stan model for counts inference.
#
# A stan file with the counts inference model is generated.
#
# @param likelihood Boolean value specifying to include the likelihood in the model, or just sampling from the priors.
# @return A string containing a stan model
# @examples
# counts_inference_code(likelihood = T)
# @export

counts_inference_code = function(likelihood = T){
  
code = "

functions {
  vector switching_process(real t,
                          vector z,
                          real[] theta //lambda_minus,lambda_plus, omega_minus, omega_plus
                          ) {
    vector[5] dzdt;
    real lambda_minus = theta[1];
    real lambda_plus = theta[2];
    // real lambda_minus = 1.2;
    // real lambda_plus = 1.2;
    real omega_minus = theta[3];
    real omega_plus = theta[4];

    dzdt[1] = lambda_minus*z[1] + omega_minus*z[2]; // mean -
    dzdt[2] = lambda_plus*z[2] + omega_plus*z[1]; // mean +
    dzdt[3] = lambda_minus*z[1] + omega_minus*z[2] + 2*lambda_minus*z[3] + 2*omega_minus*z[5]; // var -
    dzdt[4] = lambda_plus*z[2] + omega_plus*z[1] + 2*lambda_plus*z[4] + 2*omega_plus*z[5]; // var +
    dzdt[5] = (lambda_minus + lambda_plus)*z[5] + omega_plus*z[3] + omega_minus*z[4]; // covariance

    // dzdt[3] = 2*lambda_minus*z[3] + lambda_minus*z[1] + 2*omega_minus*z[5] + omega_minus*z[2];
    // dzdt[4] = 2*lambda_plus*z[4] + lambda_plus*z[2] + 2*omega_plus*z[5] + omega_plus*z[1];
    // dzdt[5] = (lambda_minus + lambda_plus)*z[5] + omega_plus*z[3] + omega_minus*z[4];

    return dzdt;
  }
}

data {
  
  int<lower=1> n_times; 
  vector[5] z0;
  real t0;
  vector[n_times] zminus;
  vector[n_times] zplus;
  array[n_times] real t;
  real <lower = 0> alpha_ln;
  real <lower = 0> beta_ln;
  real <lower = 0> alpha_lp;
  real <lower = 0> beta_lp;
  real <lower = 0> alpha_rn;
  real <lower = 0> beta_rn;
  real <lower = 0> alpha_rp;
  real <lower = 0> beta_rp;
  
}

// transformed data {
//   real x_r[0];
//   int x_i[0];
// }

parameters {
  real<lower=0> lambda_minus;
  real<lower=0> lambda_plus;
  real<lower=0> effective_switch_rate_n; // omega_plus / lambda_minus
  real<lower=0> effective_switch_rate_p; // omega_minus / lambda_plus
}

transformed parameters {
  real theta[4];
  theta[1] = lambda_minus;
  theta[2] = lambda_plus;
  theta[3] = effective_switch_rate_p*lambda_plus; // omega_minus
  theta[4] = effective_switch_rate_n*lambda_minus; // omega_plus
}

model {
  // real z_hat[n_times,5];
  array[n_times] vector[5] z_hat = ode_rk45(switching_process, z0, t0, t, theta);
  
  target += gamma_lpdf(lambda_minus | alpha_ln, beta_ln);
  target += gamma_lpdf(lambda_plus | alpha_lp, beta_lp);
  target += gamma_lpdf(effective_switch_rate_n | alpha_rn, beta_rn);
  target += gamma_lpdf(effective_switch_rate_p | alpha_rp, beta_rp);

  // z_hat = integrate_ode_rk45(switching_process, z0, t0, t, theta, x_r, x_i);

"

if(likelihood){

  code = paste0(code,  
"
  for (i in 1:n_times) {
    target += normal_lpdf(zminus[i] | z_hat[i,1], sqrt(z_hat[i,3]));
    target += normal_lpdf(zplus[i] | z_hat[i,2], sqrt(z_hat[i,4]));
  }
" )
  
}

code = paste0(code,
"
}

generated quantities {
  real pred_minus[n_times];
  real pred_plus[n_times];
  // real var_minus[n_times];
  // real var_plus[n_times];
  // real pred[n_times, 5] = integrate_ode_rk45(switching_process, z0, t0, t, theta, x_r, x_i);
  array[n_times] vector[5] pred = ode_rk45(switching_process, z0, t0, t, theta);

  for (i in 1:n_times){
    pred_minus[i] = pred[i,1];
    pred_plus[i] = pred[i,2];
    // var_minus[i] = pred[i,3];
    // var_plus[i] = pred[i,4];
  }
}

" )  
  
return(code)
  
}

# generate_code = function(max_depth){
#   
#   tree = tibble(node = "n",level = 0, delta_m = "delta_m_n", phi = "", p_switch = "", nu = "nu_n") %>% 
#     mutate(m = "m_n = nu_n*delta_m_n", sampling_nu = "nu_n ~ beta(at,bt)", sampling_phi = "")
#   
#   for(l in 1:max_depth){
#     
#     nod = tree %>% filter(level == l-1) %>% pull(node)
#     
#     new =  lapply(1:length(nod), function(i){
#       
#       new_nodes = tibble(node = c(paste0(nod[i],"n"),paste0(nod[i],"p")), level = l,
#                          rate = c("rn","rp"), 
#                          phi = c(paste0("phi_",nod[i],"n"),paste0("phi_",nod[i],"p")),
#                          sampling_phi = 
#                            c(paste0("phi_",nod[i],"n"," ~ beta(gamma,gamma)"),
#                              paste0("phi_",nod[i],"p = 1 - ","phi_",nod[i],"n"))) %>% 
#         mutate(delta_m = paste0("delta_m_",node,"=(delta_m_",nod[i],"- m_",nod[i],")*",phi), 
#                nu = paste0("nu_",node))
#       
#       if(l < max_depth){
#         
#         new_nodes =   new_nodes %>% rowwise() %>% 
#           mutate(p_switch = paste0("p_switch_",node," = ","1-exp(-",rate,"*delta_m_",node,")"),
#                  sampling_nu = paste0("nu_",
#                         node," ~ beta(p_switch_",node,"/",rate,",delta_m_",node,"- p_switch_",node,"/",rate,")")) %>% 
#           mutate(m = paste0("m_",node,"= ","nu_",node,"*","delta_m_",node))
#         
#       }else{
#         
#         new_nodes =  new_nodes %>% mutate(p_switch = "", sampling_nu = "", m = paste0("m_",node,"= delta_m_",node))   
#         
#       }
#       
#       new_nodes = new_nodes %>% dplyr::select(-rate)
#       
#     }) %>% bind_rows()
#     
#     tree = rbind(tree,new)
#     
#   }
#   
#   
#   tree = tree %>% rowwise() %>% mutate(delta_t = paste0("delta_t_",node), 
#                                        sampling_delta_t = paste0(delta_t, " ~ gamma(m_",node,",2*mu*l)"))
#   
#   max_deepth = tree$level %>% max()
#   all_nodes = tree %>% pull(node)
#   leaves = tree %>% mutate(leave = ifelse(!paste0(node,"p") %in% 
#                                             all_nodes & !paste0(node,"n") %in% all_nodes,T,F)) %>%
#     filter(leave) %>% pull(node)
#   coeff = list("1","1","(1+s)","1/(1+s)")
#   names(coeff) = c("nn","pp","np","pn")
#   
#   ccf_leaves =  lapply(leaves,function(leave){
#     
#     t = tree %>% filter(node == leave) %>% pull(delta_t)
#     last = substr(leave,start = nchar(leave),stop = nchar(leave))
#     for (j in 1:(nchar(leave)-1)){
#       ch = substr(leave,start = 1,stop = nchar(leave)-j)
#       delta_t = tree %>% filter(node == ch) %>% pull(delta_t)
#       t = paste0(t,"+",delta_t,"*",coeff[[paste0(substr(ch,start = nchar(ch),stop = nchar(ch)),last)]])
#     }
#     tibble(leave = leave,last = last, t = t, ccf = paste0("exp(-(",t,"))"))
#   }) %>% bind_rows()
#   
#   tot = ccf_leaves %>% group_by(last) %>% summarize(n = paste0(ccf,collapse = "+"))
#   ccf_minus = tot %>% filter(last == "n") %>% pull(n)
#   ccf_plus = tot %>% filter(last == "p") %>% pull(n)
#   ccf_leaves = ccf_leaves %>% mutate(ccf_minus = ifelse(last == "n",
#                                                         paste0(ccf,"/(",ccf_minus,")"),"1e-3")) %>% 
#     mutate(ccf_plus = ifelse(last == "p",paste0(ccf,"/(",ccf_plus,")"),"1e-3"))
#   
#   tree = full_join(tree,ccf_leaves %>% dplyr::rename(node = leave) %>% 
#                      dplyr::select(node,ccf_minus,ccf_plus),by = "node")
#   
#   for(lev in 1:max_deepth){
#     
#     nodes = tree %>% filter(level == max_deepth - lev,
#                             is.na(ccf_plus),
#                             is.na(ccf_minus)) %>% pull(node)
#     
#     for (nod in nodes){
#       
#       # ccf_1 = c(tree %>% filter(node == paste0(nod,"p")) %>% pull(ccf_plus),
#       #           tree %>% filter(node == paste0(nod,"p")) %>% pull(ccf_minus))
#       # ccf_2 = c(tree %>% filter(node == paste0(nod,"n")) %>% pull(ccf_plus),
#       #           tree %>% filter(node == paste0(nod,"n")) %>% pull(ccf_minus))
#       # 
#       # ccf = paste0(ccf_1,"+",ccf_2)
#       
#       ccf_1 = c(paste0("ccf_plus_",paste0(nod,"p")),paste0("ccf_minus_",paste0(nod,"p")))
#       ccf_2 = c(paste0("ccf_plus_",paste0(nod,"n")),paste0("ccf_minus_",paste0(nod,"n")))
#       ccf = paste0(ccf_1,"+",ccf_2)
#       
#       tree  = tree  %>% mutate(ccf_plus = ifelse(node == nod,ccf[1],ccf_plus),
#                                ccf_minus = ifelse(node == nod,ccf[2],ccf_minus))
#     }
#     
#   }
#   
#   tree = tree %>% mutate(ccf_minus = paste0("ccf_minus_",node,"=",ccf_minus),
#                          ccf_plus = paste0("ccf_plus_",node,"=",ccf_plus))
#   
#   
#   # data
#   
#   stick_breaking = 
#     "
# data{
# 
# real <lower=0> delta_m_n;
# real <lower=0> delta_t_n;
# real <lower=0> at;
# real<lower=0> bt;
# real <lower=0> an;
# real<lower=0> bn;
# real<lower=0> ap;
# real<lower=0> bp;
# real ms;
# real <lower=0> sigma;
# int n;
# int Nn[n];
# int Np[n];
# int DPn[n];
# int DPp[n];
# real gamma;
# real l;
# real mu;
# " 
# 
# for (i in 1:nrow(tree)){
#   
#  if(grepl(x = tree$ccf_minus[i],pattern = "1e-3")){
#   
#   stick_breaking = paste0(stick_breaking,"real <lower=0>",gsub(x = tree$ccf_minus[i],pattern = "=1e-3",""),";\n")}
# 
#  if(grepl(x = tree$ccf_plus[i],pattern = "1e-3")){ 
#   stick_breaking = paste0(stick_breaking,"real <lower=0>",gsub(x = tree$ccf_plus[i],pattern = "=1e-3",""),";\n")}
#   
# }
# 
# # params
# 
# stick_breaking = paste0(stick_breaking, 
#                         " }
#                         
#                         
# parameters{
# real <lower=0,upper = 1> rn;
# real <lower=0, upper = 1 > rp;
# real <lower=0> s;
#   "
# )
# 
# for (i in 1:nrow(tree)){
#   
#   if((i %% 2) == 0){
#     stick_breaking = paste(stick_breaking,"real <lower=0, upper=1> ",tree$phi[i],";\n")}
#   
#   if(tree$level[i] < max_depth){
#     stick_breaking = paste(stick_breaking,"real <lower=0, upper=1> ",tree$nu[i],";\n")}
#   
#   if(i > 1){
#     stick_breaking = paste(stick_breaking,"real <lower=0> ",tree$delta_t[i],";\n")}
#   
# }
# 
# # transformed params
# 
# stick_breaking = paste0(stick_breaking," } 
# 
# transformed parameters{ 
#                         ")
# 
# for (i in 1:nrow(tree)){
#   
# 
#   if((i %% 2) == 1 & i>1){
#     stick_breaking = paste(stick_breaking,"real <lower=0, upper= 1>",tree$sampling_phi[i],";\n")}
#   
#   if( i > 1){
#     
#     stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$delta_m[i],";\n")}
#   
#    stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$m[i],";\n")
#   
#  if(i > 1 & tree$level[i] < max_depth){
#     stick_breaking = paste0(stick_breaking,"real <lower=0, upper= 1>",tree$p_switch[i],";\n")}
#   
# }
# 
# for (i in 1:nrow(tree)){
#   
# k = nrow(tree) - i + 1
# 
# if(!grepl(x = tree$ccf_minus[k],pattern = "1e-3")){
# stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$ccf_minus[k],";\n")}
#  
# if(!grepl(x = tree$ccf_plus[k],pattern = "1e-3")){ 
# stick_breaking = paste0(stick_breaking,"real <lower=0>",tree$ccf_plus[k],";\n")}
#   
# }
# 
# stick_breaking = paste0(stick_breaking, "
# 
# } 
# 
# model{
#        rn ~ beta(an,bn);
#        rp ~ beta(ap,bp);
#        s ~  lognormal(ms,sigma); \n")
# 
# for (i in 1:nrow(tree)){
#   
#   if((i %% 2) == 0){
#     stick_breaking = paste(stick_breaking,tree$sampling_phi[i],";\n")}
#   
#   if(tree$level[i] < max_depth){
#     stick_breaking = paste0(stick_breaking,tree$sampling_nu[i],";\n")}
#   
#   if(i > 1){  
#     stick_breaking = paste0(stick_breaking,tree$sampling_delta_t[i],";\n")}
#   
# }
# 
# 
# stick_breaking = paste0(stick_breaking,"
# 
# for (i in 1:n){
# 
# target += log_sum_exp([log(m_n/delta_m_n) + binomial_lpmf(Nn[i] | DPn[i], 0.5*ccf_minus_n) + 
#                                             binomial_lpmf(Np[i] | DPp[i], 0.5*ccf_plus_n)\n ")
# 
# for (i in 2:nrow(tree)){
#   
#  stick_breaking = paste0(stick_breaking,"        , log(m_",tree$node[i],"/delta_m_n) + \n  ",
#                             "        binomial_lpmf(Nn[i] | DPn[i], 0.5*ccf_minus_",tree$node[i],") + \n",
#                             "        binomial_lpmf(Np[i] | DPp[i], 0.5*ccf_plus_",tree$node[i],")\n")
# 
#   
# }
# 
# 
# # generated quantities
# 
# stick_breaking = paste0(stick_breaking,
#                         " 
# ]);
#    }
# }
# 
# generated quantities{
#                         
#      real<lower=0> effective_switch_rate_n = rn*2*mu*l;
#      real<lower=0> effective_switch_rate_p = rp*2*mu*l;                   
#                         
#     }")
# 
# return(stick_breaking)
# 
# }
