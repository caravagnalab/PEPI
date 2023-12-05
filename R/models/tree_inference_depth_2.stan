
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

 }
                        
                        
parameters{
real <lower=0,upper = 1> rn;
real <lower=0,upper = 1 > rp;
   real <lower = 0,upper = 1> vaf_minus_n; 
 real <lower = 0,upper = 1> vaf_plus_n; 
 real <lower=0, upper=1>  nu_n ;
 real <lower = 0,upper = 1 > w_minus_nn; 
 real <lower = 0,upper = 1 > w_plus_nn; 
 real <lower=0, upper=1>  phi_nn ;
 real <lower=0, upper=1>  nu_nn ;
 real <lower=0, upper=1>  nu_np ;
 real <lower = 0,upper = 1> vaf_minus_nnn; 
 real <lower = 0,upper = 1> vaf_plus_nnn; 
 real <lower=0, upper=1>  phi_nnn ;
 real <lower = 0,upper = 1> vaf_minus_nnp; 
 real <lower = 0,upper = 1> vaf_plus_nnp; 
 real <lower = 0,upper = 1> vaf_minus_npn; 
 real <lower = 0,upper = 1> vaf_plus_npn; 
 real <lower=0, upper=1>  phi_npn ;
 real <lower = 0,upper = 1> vaf_minus_npp; 
 real <lower = 0,upper = 1> vaf_plus_npp; 
 } 

transformed parameters{ 
real <lower=0>m_n = nu_n*delta_m_n;
 real <lower = 0,upper = vaf_minus_n> vaf_minus_nn = w_minus_nn*vaf_minus_n; 
 real <lower = 0,upper = vaf_plus_n> vaf_plus_nn = w_plus_nn*vaf_plus_n; 
real <lower=0>delta_m_nn=(delta_m_n- m_n)*phi_nn;
real <lower=0>m_nn= nu_nn*delta_m_nn;
real <lower=0, upper= 1>p_switch_nn = 1-exp(-rn*delta_m_nn);
 real <lower = 0,upper = vaf_minus_n> vaf_minus_np = (1-w_minus_nn)*vaf_minus_n; 
 real <lower = 0,upper = vaf_plus_n> vaf_plus_np = (1-w_plus_nn)*vaf_plus_n; 
 real <lower=0, upper= 1> phi_np = 1 - phi_nn ;
real <lower=0>delta_m_np=(delta_m_n- m_n)*phi_np;
real <lower=0>m_np= nu_np*delta_m_np;
real <lower=0, upper= 1>p_switch_np = 1-exp(-rp*delta_m_np);
real <lower=0>delta_m_nnn=(delta_m_nn- m_nn)*phi_nnn;
real <lower=0>m_nnn= delta_m_nnn;
 real <lower=0, upper= 1> phi_nnp = 1 - phi_nnn ;
real <lower=0>delta_m_nnp=(delta_m_nn- m_nn)*phi_nnp;
real <lower=0>m_nnp= delta_m_nnp;
real <lower=0>delta_m_npn=(delta_m_np- m_np)*phi_npn;
real <lower=0>m_npn= delta_m_npn;
 real <lower=0, upper= 1> phi_npp = 1 - phi_npn ;
real <lower=0>delta_m_npp=(delta_m_np- m_np)*phi_npp;
real <lower=0>m_npp= delta_m_npp;


} 

model{
rn ~ beta_proportion(rate_n,qn);
rp ~ beta_proportion(rate_p,qp); 
 vaf_minus_n~ beta_proportion(0.5*rho_n,k); 
 vaf_plus_n~ beta_proportion(0.5*rho_p,k); 
nu_n ~ beta_proportion(nu_t,qt);
 w_minus_nn~ uniform(0,1); 
 w_plus_nn~ uniform(0,1); 
 phi_nn ~ beta(gamma,gamma) ;
nu_nn ~ beta(p_switch_nn/rn,delta_m_nn- p_switch_nn/rn);
nu_np ~ beta(p_switch_np/rp,delta_m_np- p_switch_np/rp);
 vaf_minus_nnn~ beta_proportion(vaf_minus_nn + 1e-3,k); 
 vaf_plus_nnn~ beta(1,1e6); 
 phi_nnn ~ beta(gamma,gamma) ;
 vaf_minus_nnp~ beta(1,1e6); 
 vaf_plus_nnp~ beta_proportion(vaf_plus_nn + 1e-3,k); 
 vaf_minus_npn~ beta_proportion(vaf_minus_np + 1e-3,k); 
 vaf_plus_npn~ beta(1,1e6); 
 phi_npn ~ beta(gamma,gamma) ;
 vaf_minus_npp~ beta(1,1e6); 
 vaf_plus_npp~ beta_proportion(vaf_plus_np + 1e-3,k); 


for (i in 1:n){

target += log_sum_exp([log(m_n/delta_m_n) + binomial_lpmf(Nn[i] | DPn[i], vaf_minus_n) + 
                                            binomial_lpmf(Np[i] | DPp[i], vaf_plus_n)
         , log(m_nn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nn)
        , log(m_np/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_np) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_np)
        , log(m_nnn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nnn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nnn)
        , log(m_nnp/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nnp) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nnp)
        , log(m_npn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_npn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_npn)
        , log(m_npp/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_npp) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_npp)
 
]);
   }

}

generated quantities{
                        
     real<lower=0> effective_switch_rate_n = rn*2*mu*l;
     real<lower=0> effective_switch_rate_p = rp*2*mu*l;                   
                        
    }
