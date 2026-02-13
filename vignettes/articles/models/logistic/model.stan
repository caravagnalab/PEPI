functions{
  
  vector Z(real lambda_n,real s,real omega_n,real omega_p,real t_1, real t_2, vector Z0){
    
    vector[2] N;
    
    real lambda_p = lambda_n*(1+s);
    real x1 = (lambda_p + lambda_n)/2 + sqrt((lambda_p - lambda_n)^2 + 4*omega_p*omega_n)/2;
    real x2 = (lambda_p + lambda_n)/2 - sqrt((lambda_p - lambda_n)^2 + 4*omega_p*omega_n)/2;
    real c1 = (2*(x1-lambda_p)*Z0[1] + 2*omega_n*Z0[2])/(4*(x1-lambda_p)^2 + 4*omega_p*omega_n);
    real c2 = (-2*(x1-lambda_p)*Z0[2] + 2*omega_p*Z0[1])/(4*(x1-lambda_p)^2 + 4*omega_p*omega_n);
    
    N[1] = 2*c1*(x1-lambda_p)*exp(x1*(t_2-t_1)) + 2*c2*omega_n*exp(x2*(t_2 - t_1));
    N[2] = 2*c1*omega_p*exp(x1*(t_2-t_1)) + 2*c2*(x2-lambda_n)*exp(x2*(t_2 - t_1));  
    
    return N;
    
  }
  
// real lpdf_waiting_time(real t,real t_start,real t_end,real omega_p,real lambda_n){
// 
//    return log(omega_p) + lambda_n*t -
//            2*log(1 + omega_p/lambda_n*exp(lambda_n*t_start)*(exp(lambda_n*(t-t_start)) - 1)) -
//           log(1 - 1/(1 + omega_p/lambda_n*exp(lambda_n*t_start)*(exp(lambda_n*(t_end-t_start)) - 1)));
// 
// }

  
  vector switching_process(real t,
  vector z,
  array[] real theta //lambda_minus_epi,lambda_plus_epi, omega_minus_epi, omega_plus
  ){
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

data{
  
  int N_clades;
  int <lower=2> n_times;
  int N_driver_n;
  int N_dc;
  int N_cd;
  
  array[N_driver_n] int <lower=0> m_driver_n;
  array[N_driver_n,n_times,2] real <lower=0> ccf_driver_n;
  array[N_driver_n] real ms_driver_n;
  array[N_driver_n] real <lower=0> sigma_driver_n;
  
  array[N_dc,2] int <lower=0> m_dc;
  array[N_dc,n_times,2] real <lower=0> ccf_dc_driver;
  array[N_dc,n_times,2] real <lower=0> ccf_dc_clade;
  array[N_dc] real ms_dc;
  array[N_dc] real <lower=0> sigma_dc;
  
  array[N_cd,2] int <lower=0> m_cd;
  array[N_cd,n_times,2] real <lower=0> ccf_cd_driver;
  array[N_cd,n_times,2] real <lower=0> ccf_cd_clade;
  array[N_cd] real ms_cd;
  array[N_cd] real <lower=0> sigma_cd;
  
  array[n_times] real <lower=0> times;
  array[n_times] real <lower=0> zminus;
  array[n_times] real <lower=0> zplus;
  array[N_clades] int <lower=0> m_clade;
  array[N_clades,n_times,2] real <lower=0> ccf_clade;
  array[n_times,2] real <lower=0> ccf_wt;
  real t_min;
  real ms_epi;
  real <lower=0> sigma_epi;
  real <lower=0> alpha_lambda;
  real <lower=0> beta_lambda;
  int <lower=0,upper=1> include_bp;
  real <lower=0> alpha_minus;
  real <lower=0> beta_minus;
  real <lower=0> alpha_plus;
  real <lower=0> beta_plus;
  real <lower=0> alpha_n;
  real <lower=0> beta_n;
  real <lower=0> alpha_p;
  real <lower=0> beta_p;
  real <lower=0> mu;
  real <lower=0> l;
  real <lower=0> ccf_thr_count;
  real <lower=0> ccf_thr_clade;
  
}

transformed data{
  
 array[N_clades,n_times,2] real <lower=0> z_clade;  
  
  for (i in 1:N_clades){
    
    for(j in 1:n_times){
      
      z_clade[i,j,1] = zminus[j]*ccf_clade[i,j,1];
      z_clade[i,j,2] = zplus[j]*ccf_clade[i,j,2];
      
    }
    
  }


array[N_dc,n_times,2] real z_dc_clade;
array[N_dc,n_times,2] real z_dc_driver;

for (i in 1:N_dc){
  
for(j in 1:n_times){

z_dc_clade[i,j,1] = zminus[j]*ccf_dc_clade[i,j,1];
z_dc_clade[i,j,2] = zplus[j]*ccf_dc_clade[i,j,2];

z_dc_driver[i,j,1] = zminus[j]*ccf_dc_driver[i,j,1];
z_dc_driver[i,j,2] = zplus[j]*ccf_dc_driver[i,j,2];

}

}




array[N_cd,n_times,2] real  z_cd_clade;
array[N_cd,n_times,2] real z_cd_driver;

for (i in 1:N_cd){
  
for(j in 1:n_times){

z_cd_clade[i,j,1] = zminus[j]*ccf_cd_clade[i,j,1];
z_cd_clade[i,j,2] = zplus[j]*ccf_cd_clade[i,j,2];

z_cd_driver[i,j,1] = zminus[j]*ccf_cd_driver[i,j,1];
z_cd_driver[i,j,2] = zplus[j]*ccf_cd_driver[i,j,2];

}

}




array[N_driver_n,n_times,2] real z_driver_n;

for (i in 1:N_driver_n){
  
for(j in 1:n_times){

z_driver_n[i,j,1] = zminus[j]*ccf_driver_n[i,j,1];
z_driver_n[i,j,2] = zplus[j]*ccf_driver_n[i,j,2];

}

}

array[n_times,2] real z_wt;

for (j in 1:n_times){

z_wt[j,1] =  zminus[j]*ccf_wt[j,1];
z_wt[j,2] =  zplus[j]*ccf_wt[j,2];

}


}

parameters{
  
  real <lower = 0> s_hat;
  real <lower = 0> sigma_minus;
  real <lower=0> sigma_plus;
  real <lower = t_min, upper=times[1]> tmrca;
  array[N_clades] real<lower=tmrca,upper=times[1]> t_clade; 
  

  array[N_driver_n] real <lower = 0> s_driver_n;
  array[N_driver_n] real<lower=tmrca,upper = times[1]> t_driver_n;

  

  array[N_dc] real <lower = 0> s_dc;
  array[N_dc,2]  real<lower=tmrca,upper = times[1]> t_dc;



  array[N_cd] real <lower = 0> s_cd;
  array[N_cd,2] real<lower=tmrca,upper = times[1]> t_cd;

  real <lower = 0> lambda_n; 
  real <lower = 0> omega_p;
  real <lower = 0> omega_n;

  
}

transformed parameters{
  
real <lower=-1> s_epi;
s_epi = s_hat - 1;

array[4] real theta;

theta[1] = lambda_n; // lambda_n
theta[2] = lambda_n*(1+s_epi); // lambda_p
theta[3] = omega_n; // omega_minus
theta[4] = omega_p; // omega_plus

// real delta = (theta[2] + theta[1])/2 + sqrt((theta[2] - theta[1])^2 + 4*theta[4]*theta[3])/2;

}


model{
  
  // prior
  
s_hat ~ lognormal(ms_epi,sigma_epi);
lambda_n ~ gamma(alpha_lambda,beta_lambda);
sigma_minus ~ gamma(alpha_minus,beta_minus);
sigma_plus ~ gamma(alpha_plus,beta_plus);
omega_n ~ gamma(alpha_n,beta_n);
omega_p ~ gamma(alpha_p,beta_p);
tmrca ~ uniform(t_min,times[1]);

if(N_driver_n > 0){
  
 for (i in 1:N_driver_n){
    
    s_driver_n[i] ~ lognormal(ms_driver_n[i],sigma_driver_n[i]);
    t_driver_n[i] ~ uniform(tmrca,times[1]);
    m_driver_n[i] ~ poisson(2*mu*l*lambda_n*2*(t_driver_n[i] - tmrca));
    

if(ccf_driver_n[i,1,1]*ccf_driver_n[i,1,2] > 0){
  
array[4] real theta_driver_n;

for (r in 1:4){
  
 theta_driver_n[r] = theta[r]*(1+s_driver_n[i]);
 
}
    
array[n_times - 1] vector[5] z_hat_n = ode_rk45(switching_process, [z_driver_n[i,1,1],z_driver_n[i,1,2],0,0,0]',
                                            times[1], times[2:n_times], theta_driver_n);

for (k in 1:(n_times-1)){

   if(ccf_driver_n[i,k+1,1] >= ccf_thr_clade){    
     target += normal_lpdf(log(z_driver_n[i,k+1,1]) | log(z_hat_n[k,1]), sigma_minus);}

   if(ccf_driver_n[i,k+1,2] >= ccf_thr_clade){   
     target += normal_lpdf(log(z_driver_n[i,k+1,2]) | log(z_hat_n[k,2]), sigma_plus);}

}

}else{
    
  for (j in 1:(n_times-1)){

  target += normal_lpdf(log(z_driver_n[i,j+1,1]) |
      log(z_driver_n[i,1,1]*exp(lambda_n*(1+s_driver_n[i])*(times[j+1] - times[1]))), sigma_minus);}
      
}

 for(j in 1:n_times){
      
      vector[2] z = Z(lambda_n*(1 + s_driver_n[i]),s_epi,omega_n*(1 + s_driver_n[i]),
                      omega_p*(1 + s_driver_n[i]),t_driver_n[i],times[j],[1,0]');

if(ccf_driver_n[i,j,1] >= ccf_thr_clade){
  
     z_driver_n[i,j,1] ~ exponential(1/z[1]);
}   
     
if(ccf_driver_n[i,j,2] > ccf_thr_clade){

      z_driver_n[i,j,2] ~ exponential(1/z[2]);
  }
      
}
    
  }
}

if(N_clades > 0){
  
   N_clades ~ poisson((omega_p/lambda_n)*z_wt[1,1]*(1/((ccf_thr_clade*z_wt[1,2])^(1/(1+s_epi))) - 1/((z_wt[1,2])^(1/(1+s_epi)))));
  
  for (i in 1:N_clades){
    
 if(i == 1){ 
              t_clade[i] ~ logistic(tmrca + log(lambda_n/omega_p)/lambda_n,1/lambda_n) T[tmrca,times[1]];
  }
   
    m_clade[i] ~ poisson(2*mu*l*lambda_n*2*(t_clade[i] - tmrca));
   
  for (j in 1:n_times){
        
        vector[2] z = Z(lambda_n,s_epi,omega_n,omega_p,t_clade[i],times[j],[0,1]');
        
      if(ccf_clade[i,j,2] > ccf_thr_clade){
        
                       z_clade[i,j,2] ~ exponential(1/z[2]);
                    
                    // z_clade[i,j,2] ~ exponential(exp(-lambda_n*(1+s_epi)*(times[j] - t_clade[i])));
         
      }
          
        if(ccf_clade[i,j,1] > ccf_thr_clade){

              z_clade[i,j,1] ~ exponential(1/z[1]);
              
           // z_clade[i,j,1] ~ exponential((delta - lambda_n)*exp(-lambda_n*(1+s_epi)*(times[j] - t_clade[i]))/omega_n);
        }
          
      }
      
      
    }
    
}
  


if(N_dc > 0){
  
  for (i in 1:N_dc){
    
 target += poisson_lpmf( 1 | (omega_p/lambda_n)*z_dc_driver[i,1,1]*(
   1/((ccf_thr_clade*z_dc_driver[i,1,2])^(1/(1+s_epi))) - 1/((z_dc_driver[i,1,2])^(1/(1+s_epi)))));
   
    s_dc[i] ~ lognormal(ms_dc[i],sigma_dc[i]);
    
    t_dc[i,1] ~ uniform(tmrca,times[1]);
    
    t_dc[i,2] ~ logistic(t_dc[i,1] + 
                        log(lambda_n/omega_p)/(lambda_n*(1 + s_dc[i])),
                         1/(lambda_n*(1 + s_dc[i]))) T[t_dc[i,1],times[1]]; 
    
    
    m_dc[i,1] ~ poisson(2*mu*l*lambda_n*2*(t_dc[i,1] - tmrca));
    m_dc[i,2] ~ poisson(2*mu*l*lambda_n*(1 + s_dc[i])*2*(t_dc[i,2] - t_dc[i,1]));
    
   if(ccf_dc_driver[i,1,1]*ccf_dc_driver[i,1,2] > 0){
        
        array[4] real theta_dc;
        
        for (r in 1:4){
          
          theta_dc[r] = theta[r]*(1+s_dc[i]);
          
        }
        
        array[n_times - 1] vector[5] z_hat = ode_rk45(switching_process, [z_dc_driver[i,1,1],z_dc_driver[i,1,2],0,0,0]',
        times[1], times[2:n_times], theta_dc);
        
        for (k in 1:(n_times-1)){
          
          if(ccf_dc_driver[i,k+1,1] > ccf_thr_clade){    
            target += normal_lpdf(log(z_dc_driver[i,k+1,1]) | log(z_hat[k,1]), sigma_minus);}
            
            if(ccf_dc_driver[i,k+1,2] > ccf_thr_clade){   
              target += normal_lpdf(log(z_dc_driver[i,k+1,2]) | log(z_hat[k,2]), sigma_plus);}
              
        }
        
      }
    
  for (j in 1:n_times){
      
      
      vector[2] z = Z(lambda_n*(1 + s_dc[i]),s_epi,omega_n*(1 + s_dc[i]),
      omega_p*(1 + s_dc[i]),t_dc[i,2],times[j],[0,1]');
      vector[2] zdr = Z(lambda_n*(1 + s_dc[i]),s_epi,omega_n*(1 + s_dc[i]),
      omega_p*(1 + s_dc[i]),t_dc[i,1],times[j],[1,0]');
      
      if(ccf_dc_clade[i,j,2]  > ccf_thr_clade){
        
        z_dc_clade[i,j,2]  ~ exponential(1/z[2]);
        
      }
        
      if(ccf_dc_clade[i,j,1]  > ccf_thr_clade){

        z_dc_clade[i,j,1]  ~ exponential(1/z[1]);
      }
      
      if( ccf_dc_driver[i,j,1] > ccf_thr_clade){
        
        z_dc_driver[i,j,1] ~ exponential(1/zdr[1]);
      }
      
       if( ccf_dc_driver[i,j,2] > ccf_thr_clade){
         
        z_dc_driver[i,j,2] ~ exponential(1/zdr[2]);
      }
      
    }
    
  }
  
}


if(N_cd > 0){
  
  N_cd ~ poisson((omega_p/lambda_n)*z_wt[1,1]*(1/((ccf_thr_clade*z_wt[1,2])^(1/(1+s_epi))) - 1/((z_wt[1,2])^(1/(1+s_epi)))));
  
  for (i in 1:N_cd){
    
    s_cd[i] ~ lognormal(ms_cd[i],sigma_cd[i]);
    
    t_cd[i,1] ~ logistic(tmrca + log(lambda_n/omega_p)/lambda_n,1/lambda_n) T[tmrca,times[1]]; 
    t_cd[i,2] ~ uniform(t_cd[i,1],times[1]);   
    
    m_cd[i,1] ~ poisson(2*mu*l*lambda_n*2*(t_cd[i,1] - tmrca));
    m_cd[i,2] ~ poisson(2*mu*l*lambda_n*(1 + s_epi)*2*(t_cd[i,2] - t_cd[i,1]));


  if(ccf_cd_driver[i,1,1]*ccf_cd_driver[i,1,2] > 0){
    
    array[4] real theta_cd;
    
    for (r in 1:4){
      
      theta_cd[r] = theta[r]*(1+s_cd[i]);
      
    }
    
    array[n_times - 1] vector[5] z_hat = ode_rk45(switching_process, [z_cd_driver[i,1,1],z_cd_driver[i,1,2],0,0,0]',
    times[1], times[2:n_times], theta_cd);
    
    for (k in 1:(n_times-1)){
      
      if(ccf_cd_driver[i,k+1,1] > ccf_thr_clade){    
        target += normal_lpdf(log(z_cd_driver[i,k+1,1]) | log(z_hat[k,1]), sigma_minus);}
        
        if(ccf_cd_driver[i,k+1,2] > ccf_thr_clade){   
          target += normal_lpdf(log(z_cd_driver[i,k+1,2]) | log(z_hat[k,2]), sigma_plus);}
          
    }
    
  }
  

    
    for (j in 1:n_times){
      
      vector[2] z = Z(lambda_n,s_epi,omega_n,omega_p,t_cd[i,1],times[j],[0,1]');
      vector[2] zdr = Z(lambda_n*(1+s_cd[i]),s_epi,omega_n*(1+s_cd[i]),omega_p*(1+s_cd[i]),t_cd[i,2],times[j],[0,1]');
      
      if(ccf_cd_clade[i,j,2]  > ccf_thr_clade){
        
        z_cd_clade[i,j,2]  ~ exponential(1/z[2]);
        
      }
      
      if(ccf_cd_clade[i,j,1]  > ccf_thr_clade){

        z_cd_clade[i,j,1]  ~ exponential(1/z[1]);
      }
      
    if(ccf_cd_driver[i,j,1] > ccf_thr_clade){
        
        z_cd_driver[i,j,1] ~ exponential(1/zdr[1]);
        
      }
      
      if(ccf_cd_driver[i,j,2] > ccf_thr_clade){
        
        z_cd_driver[i,j,2] ~ exponential(1/zdr[2]);
        
      }
      
    }
    
  }
  
}

 if(ccf_wt[1,2] > ccf_thr_count){
    
    array[n_times - 1] vector[5] z_hat = ode_rk45(switching_process, [z_wt[1,1],z_wt[1,2],0,0,0]',
    times[1], times[2:n_times], theta);

    for (i in 1:(n_times-1)){
      
      // vector[2] z_hat = Z(lambda_n,s_epi,omega_n,omega_p,times[1],times[i+1],[z_wt[1,1],z_wt[1,2]]');

      if(ccf_wt[i+1,1] > ccf_thr_count){ target += normal_lpdf(log(z_wt[i+1,1]) | log(z_hat[i,1]), sigma_minus);}

      if(ccf_wt[i+1,2] > ccf_thr_count){ target += normal_lpdf(log(z_wt[i+1,2]) | log(z_hat[i,2]), sigma_plus);}

       // if(ccf_wt[i+1,1] >= ccf_thr_count){ target += normal_lpdf(log(z_wt[i+1,1]) | log(z_hat[1]), sigma_minus);}
       // 
       // if(ccf_wt[i+1,2] >= ccf_thr_count){ target += normal_lpdf(log(z_wt[i+1,2]) | log(z_hat[2]), sigma_plus);}
       
}
    
  }else{
    
    for (j in 1:(n_times-1)){
      
      target += normal_lpdf(log(z_wt[j+1,1]) | log(z_wt[1,1]*exp(lambda_n*(times[j+1] - times[1]))), sigma_minus);
      
    }
    
    
  }
  


for (j in 1:n_times){

 vector[2] z = Z(lambda_n,s_epi,omega_n,omega_p,tmrca,times[j],[1,0]');

if(ccf_wt[j,1] > ccf_thr_count){

             target += exponential_lpdf(z_wt[j,1]| 1/z[1]);

          // target += exponential_lpdf(z_wt[j,1]| 1/(exp(lambda_n*(times[j] - tmrca)) +
          //                  omega_n/(delta -lambda_n)*z_wt[j,2]));

          // target += exponential_lpdf(z_wt[j,1]| 1/(exp(lambda_n*(times[j] - tmrca)) +
          //                  omega_n/(delta -lambda_n)*z_wt[j,2]));

           // target += exponential_lpdf(z_wt[j,1]| exp(-lambda_n*(times[j] - tmrca)));

}

if(include_bp > 0){
  
if(ccf_wt[j,2] > ccf_thr_count){

       target += exponential_lpdf(z_wt[j,2]|1/z[2]);
       
     // target += exponential_lpdf(z_wt[j,2]|(lambda_n*(s_epi)/omega_p)*exp(-lambda_n*(1+s_epi)*(times[j] - tmrca)));

}

}

}


}
  
generated quantities{

array[N_dc,n_times,2] real z_dc_clade_pred;
array[N_dc,n_times,2] real z_dc_driver_pred;
array[N_dc,2] real m_dc_pred;
array[N_dc] real t_dc_clade_pred;
array[N_dc] int N_dc_clade_pred;

for (i in 1:N_dc){
  
N_dc_clade_pred[i] = poisson_rng((omega_p/lambda_n)*z_dc_driver[i,1,1]*(1/((ccf_thr_clade*z_dc_driver[i,1,2])^(1/(1+s_epi))) - 
                                   1/((z_dc_driver[i,1,2])^(1/(1+s_epi)))));

   // t_dc_clade_pred[i] = gumbel_rng(t_dc[i,1] + log(lambda_n/omega_p)/(lambda_n*(1 + s_dc[i])),
   //                               1/(lambda_n*(1 + s_dc[i])));
    t_dc_clade_pred[i] = logistic_rng(t_dc[i,1] + log(lambda_n/omega_p)/(lambda_n*(1 + s_dc[i])),
                                 1/(lambda_n*(1 + s_dc[i])));
                                 
   m_dc_pred[i,1] = poisson_rng(2*mu*l*lambda_n*2*(t_dc[i,1] - tmrca));
   m_dc_pred[i,2] = poisson_rng(2*mu*l*lambda_n*(1 + s_dc[i])*2*(t_dc[i,2] - t_dc[i,1]));
    

for (j in 1:n_times){
      
      // vector[2] z = Z(lambda_n*(1 + s_dc[i]),s_epi,omega_n,omega_p,t_dc[i,2],times[j],[0,1]');
      // vector[2] zdr = Z(lambda_n*(1 + s_dc[i]),s_epi,omega_n,omega_p,t_dc[i,1],times[j],[1,0]');
      
      vector[2] z = Z(lambda_n*(1 + s_dc[i]),s_epi,omega_n*(1 + s_dc[i]),omega_p*(1 + s_dc[i]),t_dc[i,2],times[j],[0,1]');
      vector[2] zdr = Z(lambda_n*(1 + s_dc[i]),s_epi,omega_n*(1 + s_dc[i]),omega_p*(1 + s_dc[i]),t_dc[i,1],times[j],[1,0]');
      
  z_dc_clade_pred[i,j,2]  = exponential_rng(1/z[2]);
      
  z_dc_clade_pred[i,j,1]  = exponential_rng(1/z[1]);
      
  z_dc_driver_pred[i,j,1] = exponential_rng(1/zdr[1]);
    
  z_dc_driver_pred[i,j,2] = exponential_rng(1/zdr[2]);
      
}

}

array[N_cd,n_times,2] real z_cd_clade_pred;
array[N_cd,n_times,2] real z_cd_driver_pred;
array[N_cd,2] real m_cd_pred;
array[N_cd] real t_cd_clade_pred;

if(N_cd > 0){
int N_cd_pred;
N_cd_pred = poisson_rng((omega_p/lambda_n)*z_wt[1,1]*(1/((ccf_thr_clade*z_wt[1,2])^(1/(1+s_epi))) - 1/((z_wt[1,2])^(1/(1+s_epi)))));
}

for (i in 1:N_cd){
  
   // t_cd_clade_pred[i] = gumbel_rng(tmrca + log(lambda_n/omega_p)/lambda_n,1/lambda_n);
   t_cd_clade_pred[i] = logistic_rng(tmrca + log(lambda_n/omega_p)/lambda_n,1/lambda_n);
   m_cd_pred[i,1] = poisson_rng(2*mu*l*lambda_n*2*(t_cd[i,1] - tmrca));
   m_cd_pred[i,2] = poisson_rng(2*mu*l*lambda_n*(1 + s_epi)*2*(t_cd[i,2] - t_cd[i,1]));
    

for (j in 1:n_times){
      
      vector[2] z = Z(lambda_n,s_epi,omega_n,omega_p,t_cd[i,1],times[j],[0,1]');
      // vector[2] zdr = Z(lambda_n*(1+s_cd[i]),s_epi,omega_n,omega_p,t_cd[i,2],times[j],[0,1]');
      vector[2] zdr = Z(lambda_n*(1+s_cd[i]),s_epi,omega_n*(1+s_cd[i]),omega_p*(1+s_cd[i]),t_cd[i,2],times[j],[0,1]');
      
  z_cd_clade_pred[i,j,2]  = exponential_rng(1/z[2]);
      
  z_cd_clade_pred[i,j,1]  = exponential_rng(1/z[1]);
      
  z_cd_driver_pred[i,j,1] = exponential_rng(1/zdr[1]);
    
  z_cd_driver_pred[i,j,2] = exponential_rng(1/zdr[2]);
      
}

}

array[N_driver_n,n_times,2] real z_driver_n_pred;
array[N_driver_n] real m_driver_n_pred;

for (i in 1:N_driver_n){
  
  
  m_driver_n_pred[i] = poisson_rng(2*mu*l*lambda_n*2*(t_driver_n[i] - tmrca));
    

for (j in 1:n_times){
      
  // vector[2] zdr = Z(lambda_n*(1 + s_driver_n[i]),s_epi,omega_n,omega_p,t_driver_n[i],times[j],[1,0]');
   vector[2] zdr = Z(lambda_n*(1 + s_driver_n[i]),s_epi,omega_n*(1 + s_driver_n[i]),omega_p*(1 + s_driver_n[i]),
   t_driver_n[i],times[j],[1,0]');   
 
  z_driver_n_pred[i,j,1] = exponential_rng(1/zdr[1]);
    
  z_driver_n_pred[i,j,2] = exponential_rng(1/zdr[2]);
      
}

}

array[N_clades,n_times,2] real z_clade_pred;
array[N_clades] real m_clade_pred;
array[N_clades] real t_clade_pred;

if(N_clades > 0){
int N_clades_pred;  
N_clades_pred = poisson_rng((omega_p/lambda_n)*z_wt[1,1]*(1/((ccf_thr_clade*z_wt[1,2])^(1/(1+s_epi))) - 1/((z_wt[1,2])^(1/(1+s_epi)))));
}

for (i in 1:N_clades){
  
    // t_clade_pred[i] = gumbel_rng(tmrca + log(lambda_n/omega_p)/lambda_n,1/lambda_n);
    t_clade_pred[i] = logistic_rng(tmrca + log(lambda_n/omega_p)/lambda_n,1/lambda_n);
    m_clade_pred[i] = poisson_rng(2*mu*l*lambda_n*2*(t_clade[i] - tmrca));
    

for (j in 1:n_times){
      
  vector[2] z = Z(lambda_n,s_epi,omega_n,omega_p,t_clade[i],times[j],[0,1]');
      
 
  z_clade_pred[i,j,1] = exponential_rng(1/z[1]);
    
  z_clade_pred[i,j,2] = exponential_rng(1/z[2]);
      
}

}


array[n_times,2] real z_wild_type_pred;

for (j in 1:n_times){

   vector[2] z = Z(lambda_n,s_epi,omega_n,omega_p,tmrca,times[j],[1,0]');
   

 z_wild_type_pred[j,1] = exponential_rng(1/z[1]);
   
 z_wild_type_pred[j,2] = exponential_rng(1/z[2]);
   
}

// priors  

  real  s_epi_prior;
  real tmrca_prior;
  real <lower = 0> lambda_n_prior;
  real <lower = 0> omega_p_prior;
  real <lower = 0> omega_n_prior;
  vector[N_clades] t_clade_prior;
  
  s_epi_prior = lognormal_rng(ms_epi,sigma_epi) - 1;
  lambda_n_prior = gamma_rng(alpha_lambda,beta_lambda);
  omega_n_prior = gamma_rng(alpha_n,beta_n);
  omega_p_prior = gamma_rng(alpha_p,beta_p);
  tmrca_prior = uniform_rng(t_min,times[1]);
  

 for (i in 1:N_clades){
    
    // t_clade_prior[i] = gumbel_rng(tmrca_prior + log(lambda_n_prior/omega_p_prior)/lambda_n_prior,1/lambda_n_prior);
    t_clade_prior[i] = logistic_rng(tmrca_prior + log(lambda_n_prior/omega_p_prior)/lambda_n_prior,1/lambda_n_prior);
    
  }
  
  

  vector[N_driver_n] t_driver_n_prior;
  vector[N_driver_n] s_driver_n_prior;
  for (i in 1:N_driver_n){
    
    s_driver_n_prior[i] = lognormal_rng(ms_driver_n[i],sigma_driver_n[i]);
    t_driver_n_prior[i] = uniform_rng(tmrca_prior,times[1]);
    
    }

  vector[N_dc] s_dc_prior;
  matrix[N_dc,2] t_dc_prior;
  
  for (i in 1:N_dc){
    
    s_dc_prior[i] = lognormal_rng(ms_dc[i],sigma_dc[i]);
    
    t_dc_prior[i,1] = uniform_rng(tmrca_prior,times[1]);
    // t_dc_prior[i,2] = gumbel_rng(t_dc_prior[i,1] + 
    // log(lambda_n_prior/omega_p_prior)/(lambda_n_prior*(1 + s_dc_prior[i])),
    //    1/(lambda_n_prior*(1 + s_dc_prior[i]))); 
    t_dc_prior[i,2] = logistic_rng(t_dc_prior[i,1] + 
    log(lambda_n_prior/omega_p_prior)/(lambda_n_prior*(1 + s_dc_prior[i])),
       1/(lambda_n_prior*(1 + s_dc_prior[i]))); 
}
  

vector[N_cd] s_cd_prior;
   matrix[N_cd,2] t_cd_prior;
  for (i in 1:N_cd){
    
    s_cd_prior[i] = lognormal_rng(ms_cd[i],sigma_cd[i]);
    
    // t_cd_prior[i,1] = gumbel_rng(tmrca_prior + log(lambda_n_prior/omega_p_prior)/lambda_n_prior,1/lambda_n_prior); 
    t_cd_prior[i,1] = logistic_rng(tmrca_prior + log(lambda_n_prior/omega_p_prior)/lambda_n_prior,1/lambda_n_prior); 
    t_cd_prior[i,2] = uniform_rng(tmrca_prior,times[1]);   
    
}



}
