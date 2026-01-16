
data{

real <lower=0> delta_t_n;
real l;
real mu;
real ms;
real <lower=0> sigma;
real <lower=0> k;
real <lower=0>m_n;
real <lower=0,upper = 1>vaf_minus_n;
real <lower=0,upper = 1>vaf_plus_n;
real <lower=0>m_nn;
real <lower=0,upper = 1>vaf_minus_nn;
real <lower=0,upper = 1>vaf_plus_nn;
real <lower=0>m_np;
real <lower=0,upper = 1>vaf_minus_np;
real <lower=0,upper = 1>vaf_plus_np;
real <lower=0>m_nnn;
real <lower=0,upper = 1>vaf_minus_nnn;
real <lower=0,upper = 1>vaf_plus_nnn;
real <lower=0>m_nnp;
real <lower=0,upper = 1>vaf_minus_nnp;
real <lower=0,upper = 1>vaf_plus_nnp;
real <lower=0>m_npn;
real <lower=0,upper = 1>vaf_minus_npn;
real <lower=0,upper = 1>vaf_plus_npn;
real <lower=0>m_npp;
real <lower=0,upper = 1>vaf_minus_npp;
real <lower=0,upper = 1>vaf_plus_npp;
real <lower=0>m_nnnn;
real <lower=0,upper = 1>vaf_minus_nnnn;
real <lower=0,upper = 1>vaf_plus_nnnn;
real <lower=0>m_nnnp;
real <lower=0,upper = 1>vaf_minus_nnnp;
real <lower=0,upper = 1>vaf_plus_nnnp;
real <lower=0>m_nnpn;
real <lower=0,upper = 1>vaf_minus_nnpn;
real <lower=0,upper = 1>vaf_plus_nnpn;
real <lower=0>m_nnpp;
real <lower=0,upper = 1>vaf_minus_nnpp;
real <lower=0,upper = 1>vaf_plus_nnpp;
real <lower=0>m_npnn;
real <lower=0,upper = 1>vaf_minus_npnn;
real <lower=0,upper = 1>vaf_plus_npnn;
real <lower=0>m_npnp;
real <lower=0,upper = 1>vaf_minus_npnp;
real <lower=0,upper = 1>vaf_plus_npnp;
 }
                        
                        
parameters{

real <lower=0> s;
   real <lower=0>  delta_t_nn ;
 real <lower=0>  delta_t_np ;
 real <lower=0>  delta_t_nnn ;
 real <lower=0>  delta_t_nnp ;
 real <lower=0>  delta_t_npn ;
 real <lower=0>  delta_t_npp ;
 real <lower=0>  delta_t_nnnn ;
 real <lower=0>  delta_t_nnnp ;
 real <lower=0>  delta_t_nnpn ;
 real <lower=0>  delta_t_nnpp ;
 real <lower=0>  delta_t_npnn ;
 real <lower=0>  delta_t_npnp ;
 } 

transformed parameters{ 
                        real <lower=0>ccf_minus_npnp=2*vaf_minus_npnp;
real <lower=0>ccf_plus_npnp=exp(-(+delta_t_npn*(1+s)+delta_t_np*1+delta_t_n*(1+s)))/(exp(-(+delta_t_np*1+delta_t_n*(1+s)))+exp(-(+delta_t_nnn*(1+s)+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_nnp*1+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_npn*(1+s)+delta_t_np*1+delta_t_n*(1+s))));
real <lower=0>ccf_minus_npnn=exp(-(+delta_t_npn*1+delta_t_np*1/(1+s)+delta_t_n*1))/(exp(-(+delta_t_nnn*1+delta_t_nn*1+delta_t_n*1))+exp(-(+delta_t_nnp*1/(1+s)+delta_t_nn*1+delta_t_n*1))+exp(-(+delta_t_npn*1+delta_t_np*1/(1+s)+delta_t_n*1)));
real <lower=0>ccf_plus_npnn=2*vaf_plus_npnn;
real <lower=0>ccf_minus_nnpp=2*vaf_minus_nnpp;
real <lower=0>ccf_plus_nnpp=exp(-(+delta_t_nnp*1+delta_t_nn*(1+s)+delta_t_n*(1+s)))/(exp(-(+delta_t_np*1+delta_t_n*(1+s)))+exp(-(+delta_t_nnn*(1+s)+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_nnp*1+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_npn*(1+s)+delta_t_np*1+delta_t_n*(1+s))));
real <lower=0>ccf_minus_nnpn=exp(-(+delta_t_nnp*1/(1+s)+delta_t_nn*1+delta_t_n*1))/(exp(-(+delta_t_nnn*1+delta_t_nn*1+delta_t_n*1))+exp(-(+delta_t_nnp*1/(1+s)+delta_t_nn*1+delta_t_n*1))+exp(-(+delta_t_npn*1+delta_t_np*1/(1+s)+delta_t_n*1)));
real <lower=0>ccf_plus_nnpn=2*vaf_plus_nnpn;
real <lower=0>ccf_minus_nnnp=2*vaf_minus_nnnp;
real <lower=0>ccf_plus_nnnp=exp(-(+delta_t_nnn*(1+s)+delta_t_nn*(1+s)+delta_t_n*(1+s)))/(exp(-(+delta_t_np*1+delta_t_n*(1+s)))+exp(-(+delta_t_nnn*(1+s)+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_nnp*1+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_npn*(1+s)+delta_t_np*1+delta_t_n*(1+s))));
real <lower=0>ccf_minus_nnnn=exp(-(+delta_t_nnn*1+delta_t_nn*1+delta_t_n*1))/(exp(-(+delta_t_nnn*1+delta_t_nn*1+delta_t_n*1))+exp(-(+delta_t_nnp*1/(1+s)+delta_t_nn*1+delta_t_n*1))+exp(-(+delta_t_npn*1+delta_t_np*1/(1+s)+delta_t_n*1)));
real <lower=0>ccf_plus_nnnn=2*vaf_plus_nnnn;
real <lower=0>ccf_minus_npp=2*vaf_minus_npp;
real <lower=0>ccf_plus_npp=exp(-(+delta_t_np*1+delta_t_n*(1+s)))/(exp(-(+delta_t_np*1+delta_t_n*(1+s)))+exp(-(+delta_t_nnn*(1+s)+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_nnp*1+delta_t_nn*(1+s)+delta_t_n*(1+s)))+exp(-(+delta_t_npn*(1+s)+delta_t_np*1+delta_t_n*(1+s))));
real <lower=0>ccf_minus_npn=ccf_minus_npnp+ccf_minus_npnn;
real <lower=0>ccf_plus_npn=ccf_plus_npnp+ccf_plus_npnn;
real <lower=0>ccf_minus_nnp=ccf_minus_nnpp+ccf_minus_nnpn;
real <lower=0>ccf_plus_nnp=ccf_plus_nnpp+ccf_plus_nnpn;
real <lower=0>ccf_minus_nnn=ccf_minus_nnnp+ccf_minus_nnnn;
real <lower=0>ccf_plus_nnn=ccf_plus_nnnp+ccf_plus_nnnn;
real <lower=0>ccf_minus_np=ccf_minus_npp+ccf_minus_npn;
real <lower=0>ccf_plus_np=ccf_plus_npp+ccf_plus_npn;
real <lower=0>ccf_minus_nn=ccf_minus_nnp+ccf_minus_nnn;
real <lower=0>ccf_plus_nn=ccf_plus_nnp+ccf_plus_nnn;
real <lower=0>ccf_minus_n=ccf_minus_np+ccf_minus_nn;
real <lower=0>ccf_plus_n=ccf_plus_np+ccf_plus_nn;


} 

model{
       s ~  lognormal(ms,sigma); 
delta_t_nn ~ gamma(m_nn,2*mu*l);
delta_t_np ~ gamma(m_np,2*mu*l);
delta_t_nnn ~ gamma(m_nnn,2*mu*l);
delta_t_nnp ~ gamma(m_nnp,2*mu*l);
delta_t_npn ~ gamma(m_npn,2*mu*l);
delta_t_npp ~ gamma(m_npp,2*mu*l);
delta_t_nnnn ~ gamma(m_nnnn,2*mu*l);
delta_t_nnnp ~ gamma(m_nnnp,2*mu*l);
delta_t_nnpn ~ gamma(m_nnpn,2*mu*l);
delta_t_nnpp ~ gamma(m_nnpp,2*mu*l);
delta_t_npnn ~ gamma(m_npnn,2*mu*l);
delta_t_npnp ~ gamma(m_npnp,2*mu*l);
 vaf_minus_n ~ beta_proportion(0.5*ccf_minus_n,k);
 vaf_plus_n ~ beta_proportion(0.5*ccf_plus_n,k);
 vaf_minus_nn ~ beta_proportion(0.5*ccf_minus_nn,k);
 vaf_plus_nn ~ beta_proportion(0.5*ccf_plus_nn,k);
 vaf_minus_np ~ beta_proportion(0.5*ccf_minus_np,k);
 vaf_plus_np ~ beta_proportion(0.5*ccf_plus_np,k);
 vaf_minus_nnn ~ beta_proportion(0.5*ccf_minus_nnn,k);
 vaf_plus_nnn ~ beta_proportion(0.5*ccf_plus_nnn,k);
 vaf_minus_nnp ~ beta_proportion(0.5*ccf_minus_nnp,k);
 vaf_plus_nnp ~ beta_proportion(0.5*ccf_plus_nnp,k);
 vaf_minus_npn ~ beta_proportion(0.5*ccf_minus_npn,k);
 vaf_plus_npn ~ beta_proportion(0.5*ccf_plus_npn,k);
 vaf_plus_npp ~ beta_proportion(0.5*ccf_plus_npp,k);
 vaf_minus_nnnn ~ beta_proportion(0.5*ccf_minus_nnnn,k);
 vaf_plus_nnnp ~ beta_proportion(0.5*ccf_plus_nnnp,k);
 vaf_minus_nnpn ~ beta_proportion(0.5*ccf_minus_nnpn,k);
 vaf_plus_nnpp ~ beta_proportion(0.5*ccf_plus_nnpp,k);
 vaf_minus_npnn ~ beta_proportion(0.5*ccf_minus_npnn,k);
 vaf_plus_npnp ~ beta_proportion(0.5*ccf_plus_npnp,k);

                        
 } 
