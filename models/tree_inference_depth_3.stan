
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
 real <lower = 0,upper = 1 > w_minus_nnn; 
 real <lower = 0,upper = 1 > w_plus_nnn; 
 real <lower=0, upper=1>  phi_nnn ;
 real <lower=0, upper=1>  nu_nnn ;
 real <lower=0, upper=1>  nu_nnp ;
 real <lower = 0,upper = 1 > w_minus_npn; 
 real <lower = 0,upper = 1 > w_plus_npn; 
 real <lower=0, upper=1>  phi_npn ;
 real <lower=0, upper=1>  nu_npn ;
 real <lower=0, upper=1>  nu_npp ;
 real <lower = 0,upper = 1> vaf_minus_nnnn; 
 real <lower = 0,upper = 1> vaf_plus_nnnn; 
 real <lower=0, upper=1>  phi_nnnn ;
 real <lower = 0,upper = 1> vaf_minus_nnnp; 
 real <lower = 0,upper = 1> vaf_plus_nnnp; 
 real <lower = 0,upper = 1> vaf_minus_nnpn; 
 real <lower = 0,upper = 1> vaf_plus_nnpn; 
 real <lower=0, upper=1>  phi_nnpn ;
 real <lower = 0,upper = 1> vaf_minus_nnpp; 
 real <lower = 0,upper = 1> vaf_plus_nnpp; 
 real <lower = 0,upper = 1> vaf_minus_npnn; 
 real <lower = 0,upper = 1> vaf_plus_npnn; 
 real <lower=0, upper=1>  phi_npnn ;
 real <lower = 0,upper = 1> vaf_minus_npnp; 
 real <lower = 0,upper = 1> vaf_plus_npnp; 
 real <lower = 0,upper = 1> vaf_minus_nppn; 
 real <lower = 0,upper = 1> vaf_plus_nppn; 
 real <lower=0, upper=1>  phi_nppn ;
 real <lower = 0,upper = 1> vaf_minus_nppp; 
 real <lower = 0,upper = 1> vaf_plus_nppp; 
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
 real <lower = 0,upper = vaf_minus_nn> vaf_minus_nnn = w_minus_nnn*vaf_minus_nn; 
 real <lower = 0,upper = vaf_plus_nn> vaf_plus_nnn = w_plus_nnn*vaf_plus_nn; 
real <lower=0>delta_m_nnn=(delta_m_nn- m_nn)*phi_nnn;
real <lower=0>m_nnn= nu_nnn*delta_m_nnn;
real <lower=0, upper= 1>p_switch_nnn = 1-exp(-rn*delta_m_nnn);
 real <lower = 0,upper = vaf_minus_nn> vaf_minus_nnp = (1-w_minus_nnn)*vaf_minus_nn; 
 real <lower = 0,upper = vaf_plus_nn> vaf_plus_nnp = (1-w_plus_nnn)*vaf_plus_nn; 
 real <lower=0, upper= 1> phi_nnp = 1 - phi_nnn ;
real <lower=0>delta_m_nnp=(delta_m_nn- m_nn)*phi_nnp;
real <lower=0>m_nnp= nu_nnp*delta_m_nnp;
real <lower=0, upper= 1>p_switch_nnp = 1-exp(-rp*delta_m_nnp);
 real <lower = 0,upper = vaf_minus_np> vaf_minus_npn = w_minus_npn*vaf_minus_np; 
 real <lower = 0,upper = vaf_plus_np> vaf_plus_npn = w_plus_npn*vaf_plus_np; 
real <lower=0>delta_m_npn=(delta_m_np- m_np)*phi_npn;
real <lower=0>m_npn= nu_npn*delta_m_npn;
real <lower=0, upper= 1>p_switch_npn = 1-exp(-rn*delta_m_npn);
 real <lower = 0,upper = vaf_minus_np> vaf_minus_npp = (1-w_minus_npn)*vaf_minus_np; 
 real <lower = 0,upper = vaf_plus_np> vaf_plus_npp = (1-w_plus_npn)*vaf_plus_np; 
 real <lower=0, upper= 1> phi_npp = 1 - phi_npn ;
real <lower=0>delta_m_npp=(delta_m_np- m_np)*phi_npp;
real <lower=0>m_npp= nu_npp*delta_m_npp;
real <lower=0, upper= 1>p_switch_npp = 1-exp(-rp*delta_m_npp);
real <lower=0>delta_m_nnnn=(delta_m_nnn- m_nnn)*phi_nnnn;
real <lower=0>m_nnnn= delta_m_nnnn;
 real <lower=0, upper= 1> phi_nnnp = 1 - phi_nnnn ;
real <lower=0>delta_m_nnnp=(delta_m_nnn- m_nnn)*phi_nnnp;
real <lower=0>m_nnnp= delta_m_nnnp;
real <lower=0>delta_m_nnpn=(delta_m_nnp- m_nnp)*phi_nnpn;
real <lower=0>m_nnpn= delta_m_nnpn;
 real <lower=0, upper= 1> phi_nnpp = 1 - phi_nnpn ;
real <lower=0>delta_m_nnpp=(delta_m_nnp- m_nnp)*phi_nnpp;
real <lower=0>m_nnpp= delta_m_nnpp;
real <lower=0>delta_m_npnn=(delta_m_npn- m_npn)*phi_npnn;
real <lower=0>m_npnn= delta_m_npnn;
 real <lower=0, upper= 1> phi_npnp = 1 - phi_npnn ;
real <lower=0>delta_m_npnp=(delta_m_npn- m_npn)*phi_npnp;
real <lower=0>m_npnp= delta_m_npnp;
real <lower=0>delta_m_nppn=(delta_m_npp- m_npp)*phi_nppn;
real <lower=0>m_nppn= delta_m_nppn;
 real <lower=0, upper= 1> phi_nppp = 1 - phi_nppn ;
real <lower=0>delta_m_nppp=(delta_m_npp- m_npp)*phi_nppp;
real <lower=0>m_nppp= delta_m_nppp;


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
 w_minus_nnn~ uniform(0,1); 
 w_plus_nnn~ uniform(0,1); 
 phi_nnn ~ beta(gamma,gamma) ;
nu_nnn ~ beta(p_switch_nnn/rn,delta_m_nnn- p_switch_nnn/rn);
nu_nnp ~ beta(p_switch_nnp/rp,delta_m_nnp- p_switch_nnp/rp);
 w_minus_npn~ uniform(0,1); 
 w_plus_npn~ uniform(0,1); 
 phi_npn ~ beta(gamma,gamma) ;
nu_npn ~ beta(p_switch_npn/rn,delta_m_npn- p_switch_npn/rn);
nu_npp ~ beta(p_switch_npp/rp,delta_m_npp- p_switch_npp/rp);
 vaf_minus_nnnn~ beta_proportion(vaf_minus_nnn + 1e-6,k); 
 vaf_plus_nnnn~ beta(1,1e6); 
 phi_nnnn ~ beta(gamma,gamma) ;
 vaf_minus_nnnp~ beta(1,1e6); 
 vaf_plus_nnnp~ beta_proportion(vaf_plus_nnn + 1e-6,k); 
 vaf_minus_nnpn~ beta_proportion(vaf_minus_nnp + 1e-6,k); 
 vaf_plus_nnpn~ beta(1,1e6); 
 phi_nnpn ~ beta(gamma,gamma) ;
 vaf_minus_nnpp~ beta(1,1e6); 
 vaf_plus_nnpp~ beta_proportion(vaf_plus_nnp + 1e-6,k); 
 vaf_minus_npnn~ beta_proportion(vaf_minus_npn + 1e-6,k); 
 vaf_plus_npnn~ beta(1,1e6); 
 phi_npnn ~ beta(gamma,gamma) ;
 vaf_minus_npnp~ beta(1,1e6); 
 vaf_plus_npnp~ beta_proportion(vaf_plus_npn + 1e-6,k); 
 vaf_minus_nppn~ beta_proportion(vaf_minus_npp + 1e-6,k); 
 vaf_plus_nppn~ beta(1,1e6); 
 phi_nppn ~ beta(gamma,gamma) ;
 vaf_minus_nppp~ beta(1,1e6); 
 vaf_plus_nppp~ beta_proportion(vaf_plus_npp + 1e-6,k); 


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
        , log(m_nnnn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nnnn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nnnn)
        , log(m_nnnp/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nnnp) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nnnp)
        , log(m_nnpn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nnpn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nnpn)
        , log(m_nnpp/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nnpp) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nnpp)
        , log(m_npnn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_npnn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_npnn)
        , log(m_npnp/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_npnp) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_npnp)
        , log(m_nppn/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nppn) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nppn)
        , log(m_nppp/delta_m_n) + 
          binomial_lpmf(Nn[i] | DPn[i], vaf_minus_nppp) + 
        binomial_lpmf(Np[i] | DPp[i], vaf_plus_nppp)
 
]);
   }

}

generated quantities{
                        
     real<lower=0> effective_switch_rate_n = rn*2*mu*l;
     real<lower=0> effective_switch_rate_p = rp*2*mu*l;                   
                        
    }
