functions{
  
vector approx_Z(
    real lambda_n,
    real lambda_p,
    real omega_n,
    real omega_p,
    real t_1,
    real t_2,
    real W1,
    real chi
) {

  vector[2] Z;

  real dt = t_2 - t_1;

  real eps = 1e-10;

  // smooth spectral gap
  real delta =
    sqrt(square(lambda_p - lambda_n) + eps);

  // ratio
  real zeta =
    lambda_n / lambda_p;

  // smooth sine factor
  real sin_safe =
    sqrt(square(sin(pi() * zeta)) + eps);

  // dominant eigenvalue (exact 2x2 form)
  
    real x1 =
    0.5 * (lambda_p + lambda_n)
    + sqrt(
        0.25 * square(lambda_p - lambda_n)
        + omega_p * omega_n
    );

  // W2 (smooth log formulation)
  real log_W2 =
    (
      log(pi())
      + log(omega_p + eps)
      + log(W1 + eps)
      - log(lambda_p + eps)
      - log(sin_safe)
    ) / zeta 
    + log(chi) ;

  real W2 = exp(log_W2);

  // exponential modes
  real exp_n = exp(lambda_n * dt);
  real exp_p = exp(x1 * dt);

  // coupling term (smoothed only via delta)
  Z[1] = W1 * exp_n + (omega_n / delta) * W2 * exp_p;

  Z[2] =  W2 * exp_p;

  return Z;
  
}

vector approx_Z_clade(
    real lambda_n,
    real lambda_p,
    real omega_n,
    real omega_p,
    real t_1,
    real t_2,
    real W1
) {

  vector[2] Z;

  real dt = t_2 - t_1;

  real eps = 1e-10;

  // smooth spectral gap
  real delta =
    sqrt(square(lambda_p - lambda_n) + eps);

  // dominant eigenvalue (exact 2x2 form)
  real x1 =
    0.5 * (lambda_p + lambda_n)
    + sqrt(
        0.25 * square(lambda_p - lambda_n)
        + omega_p * omega_n
    );

  real exp_p = exp(x1 * dt);

  // coupling term (smoothed only via delta)
  Z[1] = (omega_n / delta) * W1 * exp_p;

  Z[2] =  W1 * exp_p;

  return Z;
  
}



  vector Z(real lambda_n,
           real s,
           real omega_n,
           real omega_p,
           real t1,
           real t2,
           vector Z0) {

    matrix[2,2] A;
    matrix[2,2] E;
    vector[2] z;

    real dt = t2 - t1;

    real lambda_p = lambda_n * (1 + s);

    A[1,1] = lambda_n;
    A[1,2] = omega_n;
    A[2,1] = omega_p;
    A[2,2] = lambda_p;

    E = matrix_exp(A * dt);

    z = E * Z0;

    return z;
    
  }

  
}

data{
  
  int N_clades_wt;
  int <lower=1> n_times;
  int <lower=1> n_sampling_times;
  int <lower=0> n_intermediate_times;
  array[n_sampling_times] int sampling_index;
  array[n_intermediate_times] int intermediate_index;
  int N_driver;
  int N_dc;
  int N_driver_n;
  int N_driver_p;
  int<lower=0,upper = 1> include_poisson;
  int m_trunk;
  int<lower=0,upper = 1> include_trunk;
  
  array[N_driver] int <lower=0> m_driver;
  array[N_driver,n_sampling_times,2] real <lower=0> ccf_driver;
  array[N_driver] real ms_driver;
  array[N_driver] real <lower=0> sigma_driver;
  
  array[N_driver_n] int <lower=0> m_driver_n;
  array[N_driver_n,n_sampling_times] real <lower=0> ccf_driver_n;
  array[N_driver_n] real ms_driver_n;
  array[N_driver_n] real <lower=0> sigma_driver_n;
  
  array[N_driver_p,n_sampling_times] real <lower=0> ccf_driver_p;
  array[N_driver_p] real ms_driver_p;
  array[N_driver_p] real <lower=0> sigma_driver_p;
  
  array[N_dc,2] int <lower=0> m_dc;
  array[N_dc,n_sampling_times,2] real <lower=0> ccf_dc_driver;
  array[N_dc,n_sampling_times,2] real <lower=0> ccf_dc_clade;
  array[N_dc] real ms_dc;
  array[N_dc] real <lower=0> sigma_dc;
  
  array[n_times] real <lower=0> times;
  array[n_sampling_times] real <lower=0> zminus;
  array[n_sampling_times] real <lower=0> zplus;
  array[n_intermediate_times] real <lower=0> ztot;
  array[N_clades_wt] int <lower=0> m_clade_wt;
  array[N_clades_wt,n_sampling_times,2] real <lower=0> ccf_clade_wt;
  array[n_sampling_times,2] real <lower=0> ccf_wt;
  real t_min;
  real ms_epi;
  real <lower=0> sigma_epi;
  real <lower=0> alpha_lambda;
  real <lower=0> beta_lambda;
  real <lower=0> alpha_n_wt;
  real <lower=0> beta_n_wt;
  real <lower=0> alpha_p_wt;
  real <lower=0> beta_p_wt;
  array[N_driver] real <lower=0>  alpha_n_driver;
  array[N_driver] real <lower=0>  beta_n_driver;
  array[N_driver] real <lower=0>  alpha_p_driver;
  array[N_driver] real <lower=0>  beta_p_driver;
  array[N_dc] real <lower=0>  alpha_n_dc;
  array[N_dc] real <lower=0>  beta_n_dc;
  array[N_dc] real <lower=0>  alpha_p_dc;
  array[N_dc] real <lower=0>  beta_p_dc;
  real <lower=0> mu;
  real <lower=0> l;
  // real <lower=0> alpha_kappa;
  // real <lower=0> beta_kappa;
  real <lower=0> min_kappa;
  real <lower=0> max_kappa;
  // real <lower=0> alpha_sigma_count;
  // real <lower=0> beta_sigma_count;
   real <lower=0> max_sigma_count;
   real <lower=0> min_sigma_count;
  
}

transformed data{
  
  real tmax = times[1];
  
  array[N_clades_wt,n_sampling_times,2] real z_clade_wt;

  array[N_dc,n_sampling_times,2] real z_dc_clade;
  array[N_dc,n_sampling_times,2] real z_dc_driver;

  array[N_driver,n_sampling_times,2] real z_driver;
  
  array[N_driver_n,n_sampling_times] real z_driver_n;
  
  array[N_driver_p,n_sampling_times] real z_driver_p;

  array[n_sampling_times,2] real z_wt;

  // frazioni
  array[n_sampling_times] real frac_wt_neg;
  array[n_sampling_times] real frac_wt_pos;

  array[N_clades_wt, n_sampling_times] real frac_clade_wt_neg;
  array[N_clades_wt, n_sampling_times] real frac_clade_wt_pos;

  array[N_dc, n_sampling_times] real frac_dc_clade_neg;
  array[N_dc, n_sampling_times] real frac_dc_clade_pos;

  array[N_dc, n_sampling_times] real frac_dc_driver_neg;
  array[N_dc, n_sampling_times] real frac_dc_driver_pos;

  array[N_driver, n_sampling_times] real frac_driver_neg;
  array[N_driver, n_sampling_times] real frac_driver_pos;


  for (i in 1:N_clades_wt){
    for(j in 1:n_sampling_times){

      z_clade_wt[i,j,1] = zminus[j] * ccf_clade_wt[i,j,1];
      z_clade_wt[i,j,2] = zplus[j]  * ccf_clade_wt[i,j,2];

    }
  }

  for (i in 1:N_dc){
    for(j in 1:n_sampling_times){

      z_dc_clade[i,j,1] = zminus[j] * ccf_dc_clade[i,j,1];
      z_dc_clade[i,j,2] = zplus[j]  * ccf_dc_clade[i,j,2];

      z_dc_driver[i,j,1] = zminus[j] * ccf_dc_driver[i,j,1];
      z_dc_driver[i,j,2] = zplus[j]  * ccf_dc_driver[i,j,2];

    }
  }

  for (i in 1:N_driver){
    for(j in 1:n_sampling_times){

      z_driver[i,j,1] = zminus[j] * ccf_driver[i,j,1];
      z_driver[i,j,2] = zplus[j]  * ccf_driver[i,j,2];

    }
  }
  
  for (i in 1:N_driver_n){
    
    for(j in 1:n_sampling_times){
      
      z_driver_n[i,j] = zminus[j] * ccf_driver_n[i,j];
      
    }
    
  }
  
    for (i in 1:N_driver_p){
    
    for(j in 1:n_sampling_times){
      
      z_driver_p[i,j] = zplus[j] * ccf_driver_p[i,j];
      
    }
    
  }

  for (j in 1:n_sampling_times){

    z_wt[j,1] = zminus[j] * ccf_wt[j,1];
    z_wt[j,2] = zplus[j]  * ccf_wt[j,2];

    frac_wt_neg[j] = z_wt[j,1] / (z_wt[j,1] + z_wt[j,2]);
    frac_wt_pos[j] = z_wt[j,2] / (z_wt[j,1] + z_wt[j,2]);

    for (i in 1:N_clades_wt){

      frac_clade_wt_neg[i,j] =
        z_clade_wt[i,j,1] / (z_clade_wt[i,j,1] + z_clade_wt[i,j,2] );

      frac_clade_wt_pos[i,j] =
        z_clade_wt[i,j,2] / (z_clade_wt[i,j,1] + z_clade_wt[i,j,2] );
    }

    for (i in 1:N_dc){

      frac_dc_clade_neg[i,j] =
        z_dc_clade[i,j,1] / (z_dc_clade[i,j,1] + z_dc_clade[i,j,2] );

      frac_dc_clade_pos[i,j] =
        z_dc_clade[i,j,2] / (z_dc_clade[i,j,1] + z_dc_clade[i,j,2] );

      frac_dc_driver_neg[i,j] =
        z_dc_driver[i,j,1] / (z_dc_driver[i,j,1] + z_dc_driver[i,j,2] );

      frac_dc_driver_pos[i,j] =
        z_dc_driver[i,j,2] / (z_dc_driver[i,j,1] + z_dc_driver[i,j,2] );
    }

    for (i in 1:N_driver){

      frac_driver_neg[i,j] =
        z_driver[i,j,1] / (z_driver[i,j,1] + z_driver[i,j,2]);

      frac_driver_pos[i,j] =
        z_driver[i,j,2] / (z_driver[i,j,1] + z_driver[i,j,2]);
    }
  }

}

parameters{
  
  real <lower = 0> lambda_n;
  real <lower = 0> s_epi;
  real <lower = 0> omega_p_wt;
  real <lower = 0> omega_n_wt;
  real <lower = 0> bc_wt;
  real<lower=0, upper=pi()> U_wt;  
  real<lower=0> W_wt;
  
  
  real <lower = t_min, upper=tmax> tmrca;
  array[N_clades_wt] real<lower=tmrca,upper=times[sampling_index[1]]> t_clade_wt; 
  array[N_clades_wt] real<lower=0> bc_clade_wt;
  
  array[N_driver] real <lower = 0> s_driver;
  array[N_driver] real<lower=tmrca,upper = tmax> t_driver;
  array[N_driver] real<lower=0> bc_driver;
  array[N_driver] real <lower = 0,upper=pi()> U_driver;
  array[N_driver] real <lower = 0> W_driver;
  
  array[N_driver_n] real <lower = 0> s_driver_n;
  array[N_driver_n] real<lower=tmrca,upper = tmax> t_driver_n;
  array[N_driver_n] real<lower=0> bc_driver_n;
  
  array[N_driver_p] real <lower = 0> s_driver_p;
  array[N_driver_p] real<lower=tmrca,upper = tmax> t_driver_p;
  array[N_driver_p] real<lower=0> bc_driver_p;

  array[N_dc] real <lower = 0> s_dc;
  array[N_dc,2] real <lower = 0> delta_t_dc;
  array[N_dc] real<lower=0> bc_driver_dc;
  array[N_dc] real<lower=0> bc_clade_dc;
  array[N_dc] real <lower = 0,upper=pi()> U_dc;
  array[N_dc] real <lower = 0> W_dc;

  array[N_dc] real <lower = 0> omega_p_dc;
  array[N_dc] real <lower = 0> omega_n_dc;
  array[N_driver] real <lower = 0> omega_p_driver;
  array[N_driver] real <lower = 0> omega_n_driver;
  
  real <lower = min_kappa, upper = max_kappa> kappa;
  real <lower = min_sigma_count, upper = max_sigma_count> sigma_count;


}

transformed parameters {

real alpha = 1/(1+s_epi + 1e-6);

real A_wt = (sin((1.0 - alpha) * U_wt) * pow(sin(alpha * U_wt), 
              alpha / (1.0 - alpha))) / pow(sin(U_wt), 1.0 / (1.0 - alpha));

real chi_wt = pow(A_wt / W_wt, (1.0 - alpha) / alpha);


array[N_driver] real A_driver;
array[N_driver] real chi_driver;

for (j in 1:N_driver){
  
  A_driver[j] = (sin((1.0 - alpha) * U_driver[j]) * pow(sin(alpha * U_driver[j]), 
  alpha / (1.0 - alpha))) / pow(sin(U_driver[j]), 1.0 / (1.0 - alpha));
  
  chi_driver[j] = pow(A_driver[j]/ W_driver[j], (1.0 - alpha) / alpha);
  
}


array[N_dc] real A_dc;
array[N_dc] real chi_dc;

for (j in 1:N_dc){
  
  A_dc[j] = (sin((1.0 - alpha) * U_dc[j]) * pow(sin(alpha * U_dc[j]), 
  alpha / (1.0 - alpha))) / pow(sin(U_dc[j]), 1.0 / (1.0 - alpha));
  
  chi_dc[j] = pow(A_dc[j]/ W_dc[j], (1.0 - alpha) / alpha);
  
}
  
  
  array[N_dc] ordered[2] t_dc;
  
  for (i in 1:N_dc){
    t_dc[i][1] = tmrca + delta_t_dc[i,1]; 
    t_dc[i][2] = t_dc[i][1] + delta_t_dc[i,2];
  }
  
  array[n_times] vector[2] Z_wt_lat;
  array[N_driver, n_times] vector[2] Z_driver_lat;
  array[N_driver_n, n_times] real log_Z_driver_n_lat;
  array[N_driver_p, n_times] real log_Z_driver_p_lat;
  array[N_dc, n_times] vector[2] Z_driver_dc_lat;
  array[N_dc, n_sampling_times] vector[2] Z_clade_dc_lat;
  array[N_clades_wt, n_sampling_times] vector[2] Z_clade_wt_lat;
  
  for (t in 1:n_times) {
    // WT
    // if (t == 1) {
      Z_wt_lat[t] = approx_Z(lambda_n, lambda_n*(1+s_epi), omega_n_wt, omega_p_wt,
                               tmrca, times[t], bc_wt, chi_wt);
    // } else {
    //   Z_wt_lat[t] = Z(lambda_n, s_epi, omega_n_wt, omega_p_wt, times[t-1], times[t], Z_wt_lat[t-1]);
    // }
    
    // DRIVER
    for (i in 1:N_driver) {
      // if (t == 1) {
        Z_driver_lat[i,t] = approx_Z(lambda_n*(1 + s_driver[i]),lambda_n*(1 + s_driver[i])*(1 + s_epi),
                            omega_n_driver[i], omega_p_driver[i], t_driver[i], 
                            times[t], bc_driver[i], chi_driver[i]);
      // } else {
      //   Z_driver_lat[i,t] = Z(lambda_n*(1 + s_driver[i]), s_epi, omega_n_driver[i], omega_p_driver[i],
      //        times[t-1], times[t], Z_driver_lat[i,t-1]);
      // }
    }
    
    for (i in 1:N_driver_n) {
      
       if (t == 1) {
        log_Z_driver_n_lat[i,t] =
        lambda_n * (1 + s_driver_n[i]) * (times[t] - t_driver_n[i]) 
        + log(bc_driver_n[i]);
      } else {
        log_Z_driver_n_lat[i,t] =
        log_Z_driver_n_lat[i,t-1] +
        lambda_n * (1 + s_driver_n[i]) * (times[t] - times[t-1]);
      }
      
    }
    
       for (i in 1:N_driver_p) {
      
       if (t == 1) {
        log_Z_driver_p_lat[i,t] =
        lambda_n * (1 + s_driver_p[i]) * (1 + s_epi)*(times[t] - t_driver_p[i]) 
        + log(bc_driver_p[i]);
      } else {
        log_Z_driver_p_lat[i,t] =
        log_Z_driver_p_lat[i,t-1] +
        lambda_n * (1 + s_driver_p[i]) * (1 + s_epi)*(times[t] - times[t-1]);
      }
      
    }
    
    // DC
    for (i in 1:N_dc) {
       // if (t == 1) {
        Z_driver_dc_lat[i,t] = approx_Z(lambda_n*(1 + s_dc[i]), 
                                lambda_n*(1 + s_dc[i])*(1+s_epi), 
                                omega_n_dc[i], omega_p_dc[i],
                                t_dc[i][1], times[t], bc_driver_dc[i], chi_dc[i]
                                );
      // } else {
      //   Z_driver_dc_lat[i,t] = Z(lambda_n*(1 + s_dc[i]), s_epi, omega_n_dc[i], omega_p_dc[i],
      //   times[t-1], times[t], Z_driver_dc_lat[i,t-1]);
      // }
    }
    
 
  }
  
    for (t in 1:n_sampling_times) {
      
      for (i in 1:N_dc) {
        
        
        // Z_clade_dc_lat[i,t] =  approx_Z_clade(lambda_n*(1 + s_dc[i]), 
        //           lambda_n*(1 + s_dc[i])*(1 + s_epi),
        //           omega_n_dc[i], omega_p_dc[i], t_dc[i][2], 
        //           times[t], bc_clade_dc[i]);
        
      if (t == 1) {

        Z_clade_dc_lat[i,t] = Z(lambda_n*(1 + s_dc[i]),s_epi,
                                omega_n_dc[i], omega_p_dc[i],
                                t_dc[i][2], times[sampling_index[t]], [0,bc_clade_dc[i]]');
      } else {

        Z_clade_dc_lat[i,t] = Z(lambda_n*(1 + s_dc[i]), s_epi, omega_n_dc[i], omega_p_dc[i],
        times[sampling_index[t-1]], times[sampling_index[t]], Z_clade_dc_lat[i,t-1]);

      }
  }
    
       // CLADES WT
    for (i in 1:N_clades_wt) {
      
      // Z_clade_wt_lat[i,t] =  approx_Z_clade(lambda_n, lambda_n*(1 + s_epi),
      //                        omega_n_wt, omega_p_wt, t_clade_wt[i], 
      //                        times[t], bc_clade_wt[i]);
                            
       if (t == 1) {
         Z_clade_wt_lat[i,t] =  Z(lambda_n, s_epi, omega_n_wt, omega_p_wt,
                                 t_clade_wt[i], times[sampling_index[t]], [0,bc_clade_wt[i]]');


      } else {
        Z_clade_wt_lat[i,t] = Z(lambda_n, s_epi, omega_n_wt, omega_p_wt,
        times[sampling_index[t-1]], times[sampling_index[t]], Z_clade_wt_lat[i,t-1]);
      }
    }
      
 }

  array[n_times] vector[2] frac_wt_lat;
  array[N_driver, n_times] vector[2] frac_driver_lat;
  array[N_dc, n_times] vector[2] frac_driver_dc_lat;
  array[N_clades_wt, n_sampling_times] vector[2] frac_clade_wt_lat;
  array[N_dc, n_sampling_times] vector[2] frac_clade_dc_lat;

  // Stabilizzazione e calcolo frazioni
  for (t in 1:n_times) {
    frac_wt_lat[t] = Z_wt_lat[t] / (sum(Z_wt_lat[t]) + 1e-20);
    frac_wt_lat[t] = frac_wt_lat[t] * (1 - 2*1e-10) + 1e-10;

    for (i in 1:N_driver) {
      frac_driver_lat[i,t] = Z_driver_lat[i,t] / (sum(Z_driver_lat[i,t]) + 1e-20);
      frac_driver_lat[i,t] = frac_driver_lat[i,t] * (1 - 2*1e-10) + 1e-10;
    }

    for (i in 1:N_dc) {
      frac_driver_dc_lat[i,t] = Z_driver_dc_lat[i,t] / (sum(Z_driver_dc_lat[i,t]) + 1e-20);
      frac_driver_dc_lat[i,t] = frac_driver_dc_lat[i,t] * (1 - 2*1e-10) + 1e-10;
    }

  }
  
  for (t in 1:n_sampling_times) {

    for (i in 1:N_dc) {
      
      frac_clade_dc_lat[i,t] = Z_clade_dc_lat[i,t] / (sum(Z_clade_dc_lat[i,t]) + 1e-20);
      frac_clade_dc_lat[i,t] = frac_clade_dc_lat[i,t] * (1 - 2*1e-10) + 1e-10;
    }

    for (i in 1:N_clades_wt) {
      frac_clade_wt_lat[i,t] = Z_clade_wt_lat[i,t] / (sum(Z_clade_wt_lat[i,t]) + 1e-20);
      frac_clade_wt_lat[i,t] = frac_clade_wt_lat[i,t] * (1 - 2*1e-10) + 1e-10;
    }
  }

array[n_intermediate_times] real ztot_lat;


if (n_intermediate_times > 0) {


for (k in 1:n_intermediate_times) {

  real total_neg = 0;
  real total_pos = 0;

  total_neg += Z_wt_lat[intermediate_index[k]][1];
  total_pos += Z_wt_lat[intermediate_index[k]][2];

  for (i in 1:N_driver) {
    total_neg += Z_driver_lat[i,intermediate_index[k]][1];
    total_pos += Z_driver_lat[i,intermediate_index[k]][2];
  }

  for (i in 1:N_driver_n) {
    total_neg += exp(log_Z_driver_n_lat[i,intermediate_index[k]]);
  }

  for (i in 1:N_dc) {
    total_neg += Z_driver_dc_lat[i,intermediate_index[k]][1];
    total_pos += Z_driver_dc_lat[i,intermediate_index[k]][2];
  }

  ztot_lat[k] = total_neg + total_pos;
  
}

}

}


model{

  
// kappa ~ gamma(alpha_kappa, beta_kappa);

kappa ~ uniform(min_kappa, max_kappa);

// sigma_count ~ gamma(alpha_sigma_count, beta_sigma_count);

sigma_count ~  uniform(min_sigma_count, max_sigma_count);

tmrca ~ uniform(t_min, tmax);

// wt
s_epi ~ lognormal(ms_epi, sigma_epi);
lambda_n ~ gamma(alpha_lambda, beta_lambda);
bc_wt ~ exponential(1);
U_wt ~ uniform(0.0, pi());
W_wt ~ exponential(1.0);

omega_n_wt ~ gamma(alpha_n_wt, beta_n_wt);
omega_p_wt ~ gamma(alpha_p_wt, beta_p_wt);


if(include_trunk > 0){
  
target += poisson_lccdf(m_trunk | 2*mu*l*lambda_n*2*(tmrca-t_min));

}


if(n_intermediate_times > 0){
  
for (k in 1:n_intermediate_times){

 target += normal_lpdf(log(ztot[k]) | log(ztot_lat[k]), sigma_count );

}

}

for (t in 1:n_sampling_times) {

  frac_wt_pos[t] ~ beta_proportion(
    frac_wt_lat[sampling_index[t]][2], kappa
  );
  
    target += normal_lpdf(log(z_wt[t,1] + z_wt[t,2]) | 
    log(Z_wt_lat[sampling_index[t],1] + Z_wt_lat[sampling_index[t],2] ), 
    sigma_count );

}

if (N_clades_wt > 0) {

  for (i in 1:N_clades_wt) {

     if(i == 1){ 
              t_clade_wt[i] ~ logistic(tmrca + log(lambda_n/omega_p_wt)/lambda_n,1/lambda_n) 
              T[tmrca,times[sampling_index[1]]];
  }else{
    
    t_clade_wt[i] ~ uniform(tmrca,times[sampling_index[1]]);
  }
   
    m_clade_wt[i] ~ poisson(2*mu*l*lambda_n*2*(t_clade_wt[i] - tmrca));
    
    bc_clade_wt[i] ~ exponential(1);

for (t in 1:n_sampling_times) {

 if(z_clade_wt[i,t,1] > 0){

   frac_clade_wt_pos[i,t] ~ beta_proportion(
     frac_clade_wt_lat[i,t][2], kappa
     );

     target += normal_lpdf(log(z_clade_wt[i,t,1] + z_clade_wt[i,t,2]) |
     log(Z_clade_wt_lat[i,t][1] + Z_clade_wt_lat[i,t][2] ), sigma_count );

 }else{

   target += normal_lpdf(log(z_clade_wt[i,t,2]) | log(Z_clade_wt_lat[i,t][2]), sigma_count );

 }

    }
    
  }
  
}

if (N_driver > 0) {
  
  for (i in 1:N_driver) {
    
    s_driver[i] ~ lognormal(ms_driver[i], sigma_driver[i]);
    omega_n_driver[i] ~ gamma(alpha_n_driver[i], beta_n_driver[i]);
    omega_p_driver[i] ~ gamma(alpha_p_driver[i], beta_p_driver[i]);
    
    t_driver[i] ~ uniform(tmrca, tmax);
    if(include_poisson == 1){
      
      m_driver[i] ~ poisson(4 * mu * l * lambda_n * (t_driver[i] - tmrca));
    }
    
    bc_driver[i] ~ exponential(1);
    
    U_driver[i] ~ uniform(0.0, pi());
    W_driver[i] ~ exponential(1.0);
    
    
    for (t in 1:n_sampling_times) {
      
        frac_driver_neg[i,t] ~ beta_proportion(frac_driver_lat[i,sampling_index[t]][1], kappa);
        
        target += normal_lpdf(log(z_driver[i,t,1] + z_driver[i,t,2] ) | 
         log(Z_driver_lat[i,sampling_index[t]][1] + Z_driver_lat[i,sampling_index[t]][2] ), sigma_count );


    }
  }
  
}


if (N_driver_p > 0) {
  
  for (i in 1:N_driver_p) {
    
    s_driver_p[i] ~ lognormal(ms_driver_p[i], sigma_driver_p[i]);
    
    t_driver_p[i] ~ uniform(tmrca, tmax);
    
    bc_driver_p[i] ~ exponential(1);
    
    for (t in 1:n_sampling_times) {
      
      target += normal_lpdf(log(z_driver_p[i,t]) | log_Z_driver_p_lat[i,sampling_index[t]], sigma_count );

    }
  }
  
}

if (N_driver_n > 0) {
  
  for (i in 1:N_driver_n) {
    
    s_driver_n[i] ~ lognormal(ms_driver_n[i], sigma_driver_n[i]);
    
    t_driver_n[i] ~ uniform(tmrca, tmax);
    
    if(include_poisson == 1){
      
       m_driver_n[i] ~ poisson(4 * mu * l * lambda_n * (t_driver_n[i] - tmrca));
    
    }
    
    bc_driver_n[i] ~ exponential(1);
    
    for (t in 1:n_sampling_times) {
      
      target += normal_lpdf(log(z_driver_n[i,t]) | 
      log_Z_driver_n_lat[i,sampling_index[t]], sigma_count );

    }
  }
  
}


if (N_dc > 0) {

  for (i in 1:N_dc) {

    s_dc[i] ~ lognormal(ms_dc[i], sigma_dc[i]);
    omega_n_dc[i] ~ gamma(alpha_n_dc[i], beta_n_dc[i]);
    omega_p_dc[i] ~ gamma(alpha_p_dc[i], beta_p_dc[i]);
    
    delta_t_dc[i,1] ~ uniform(0,tmax - tmrca);

    delta_t_dc[i,2] ~ logistic(
                        log(lambda_n*(1 + s_dc[i])/ omega_p_dc[i])/(lambda_n*(1 + s_dc[i])),
                         1/(lambda_n*(1 + s_dc[i]))) T[0,tmax - t_dc[i][1]]; 
    
    
    m_dc[i,1] ~ poisson(2*mu*l*lambda_n*2*(t_dc[i][1] - tmrca));
    m_dc[i,2] ~ poisson(2*mu*l*lambda_n*(1 + s_dc[i])*2*(t_dc[i][2] - t_dc[i][1]));
    
    bc_clade_dc[i] ~ exponential(1);
    bc_driver_dc[i] ~ exponential(1);
    
    U_dc[i] ~ uniform(0.0, pi());
    W_dc[i] ~ exponential(1.0);
    

    for (t in 1:n_sampling_times) {

      frac_dc_driver_pos[i,t] ~ beta_proportion(
        frac_driver_dc_lat[i,sampling_index[t]][2], kappa
      );

   target += normal_lpdf(log(z_dc_driver[i,t,1] + z_dc_driver[i,t,2] ) | 
   log(Z_driver_dc_lat[i,sampling_index[t]][1] + Z_driver_dc_lat[i,sampling_index[t]][2]), sigma_count );

 if(z_dc_clade[i,t,1] > 0){

     frac_dc_clade_pos[i,t] ~ beta_proportion(
        frac_clade_dc_lat[i,t][2], kappa
      );

      target += normal_lpdf(log(z_dc_clade[i,t,1] + z_dc_clade[i,t,2]) | 
      log(Z_clade_dc_lat[i,t][1] + Z_clade_dc_lat[i,t][2] ), sigma_count );

 }else{

   target += normal_lpdf(log(z_dc_clade[i,t,2]) | log(Z_clade_dc_lat[i,t][2]), sigma_count );

 }

    }

    }
  }
  

}
  
generated quantities {
    
    // =========================================================
    // PRIOR
    // =========================================================
 
    real lambda_n_prior = gamma_rng(alpha_lambda, beta_lambda);
    real tmrca_prior = uniform_rng(t_min, tmax);
    
    real s_epi_prior = lognormal_rng(ms_epi, sigma_epi);
    
    real omega_n_wt_prior = gamma_rng(alpha_n_wt, beta_n_wt);
    real omega_p_wt_prior = gamma_rng(alpha_p_wt, beta_p_wt);
    
    vector[N_driver] s_driver_prior;
    vector[N_driver] omega_n_driver_prior;
    vector[N_driver] omega_p_driver_prior;
    vector[N_driver_n] s_driver_n_prior;
    vector[N_driver_p] s_driver_p_prior;
    
    vector[N_dc] s_dc_prior;
    vector[N_dc] omega_n_dc_prior;
    vector[N_dc] omega_p_dc_prior;
    
    
    for (i in 1:N_driver) {
      
      s_driver_prior[i] =
      lognormal_rng(ms_driver[i], sigma_driver[i]);

      omega_n_driver_prior[i] =
      gamma_rng(alpha_n_driver[i], beta_n_driver[i]);
      
      omega_p_driver_prior[i] =
      gamma_rng(alpha_p_driver[i], beta_p_driver[i]);
    }
    
    for (i in 1:N_driver_n) {
      
      s_driver_n_prior[i] =
      lognormal_rng(ms_driver_n[i], sigma_driver_n[i]);
      
  }
  
    for (i in 1:N_driver_p) {
      
      s_driver_p_prior[i] =
      lognormal_rng(ms_driver_p[i], sigma_driver_p[i]);
      
  }
    
    
    for (i in 1:N_dc) {
      
      s_dc_prior[i] =
      lognormal_rng(ms_dc[i], sigma_dc[i]);

      omega_n_dc_prior[i] =
      gamma_rng(alpha_n_dc[i], beta_n_dc[i]);
      
      omega_p_dc_prior[i] =
      gamma_rng(alpha_p_dc[i], beta_p_dc[i]);
    }
    
    
    // =========================================================
    // TIMES PRIOR
    // =========================================================
    
    array[N_dc,2] real t_dc_prior;
    vector[N_driver] t_driver_prior;
    vector[N_driver_n] t_driver_n_prior;
    vector[N_driver_p] t_driver_p_prior;
    vector[N_clades_wt] t_clade_wt_prior;
    
    for (i in 1:N_clades_wt) {
      
      t_clade_wt_prior[i] =
      logistic_rng(
        tmrca_prior + log(lambda_n_prior / omega_p_wt_prior) / lambda_n_prior,
        1 / lambda_n_prior
        );
    }
    
    for (i in 1:N_driver) {
      
      t_driver_prior[i] =
      uniform_rng(tmrca_prior, tmax);
    }
    
     for (i in 1:N_driver_n) {
      
      t_driver_n_prior[i] =
      uniform_rng(tmrca_prior, tmax);
    }
    
     for (i in 1:N_driver_p) {
      
      t_driver_p_prior[i] =
      uniform_rng(tmrca_prior, tmax);
    }
    
    for (i in 1:N_dc) {
      
      t_dc_prior[i,1] =
      uniform_rng(tmrca_prior, tmax);
      
      t_dc_prior[i,2] =
      logistic_rng(
        t_dc_prior[i,1]
        + log(lambda_n_prior * (1 + s_dc_prior[i]) / omega_p_dc_prior[i])
        / (lambda_n_prior * (1 + s_dc_prior[i])),
        1 / (lambda_n_prior * (1 + s_dc_prior[i]))
        );
        
    }
  
    // =========================================================
    // Z POSTERIOR PREDICTIVE (t = 1)
    // =========================================================

    array[n_sampling_times,2] real z_wt_pred;
    array[N_driver,n_sampling_times,2] real z_driver_pred;
    array[N_dc,n_sampling_times,2] real z_dc_driver_pred;
    array[N_driver_n,n_sampling_times] real z_driver_n_pred;
    array[N_driver_p,n_sampling_times] real z_driver_p_pred;
    array[N_dc,n_sampling_times,2] real z_dc_clade_pred;
    array[N_clades_wt,n_sampling_times,2] real z_clade_wt_pred;

for (k in 1:n_sampling_times){

    z_wt_pred[k,1] = exp(normal_rng(log(Z_wt_lat[sampling_index[k],1]), sigma_count));
    z_wt_pred[k,2] = exp(normal_rng(log(Z_wt_lat[sampling_index[k],2]),sigma_count));

for (i in 1:N_driver_n) {
  
  
  z_driver_n_pred[i,k] = exp(normal_rng(log_Z_driver_n_lat[i,sampling_index[k]],sigma_count));
  
}

for (i in 1:N_driver_p) {
  
  
  z_driver_p_pred[i,k] = exp(normal_rng(log_Z_driver_p_lat[i,sampling_index[k]],sigma_count));
  
}

    for (i in 1:N_dc) {

      z_dc_driver_pred[i,k,1] =
      exp(normal_rng(log(Z_driver_dc_lat[i,sampling_index[k]][1]),sigma_count));

      z_dc_driver_pred[i,k,2] =
      exp(normal_rng(log(Z_driver_dc_lat[i,sampling_index[k]][2]),sigma_count));

      z_dc_clade_pred[i,k,1] =
      exp(normal_rng(log(Z_clade_dc_lat[i,k][1]),sigma_count));

      z_dc_clade_pred[i,k,2] =
      exp(normal_rng(log(Z_clade_dc_lat[i,k][2]),sigma_count));
    }

    for (i in 1:N_clades_wt){

      z_clade_wt_pred[i,k,1] =
      exp(normal_rng(log(Z_clade_wt_lat[i,k][1]),sigma_count));

      z_clade_wt_pred[i,k,2] =
      exp(normal_rng(log(Z_clade_wt_lat[i,k][2]),sigma_count));
    }
    
     for (i in 1:N_driver){

      z_driver_pred[i,k,1] =
       exp(normal_rng(log(Z_driver_lat[i,sampling_index[k]][1]),sigma_count));
      
        z_driver_pred[i,k,2] =
       exp(normal_rng(log(Z_driver_lat[i,sampling_index[k]][2]),sigma_count));
       
    }
  
}    
    // ==============================
    // FRACTIONS POSTERIOR PREDICTIVE
    // ==============================

    array[n_sampling_times] simplex[2] frac_wt_pred;
    array[N_driver, n_sampling_times] simplex[2] frac_driver_pred;
    array[N_dc, n_sampling_times] simplex[2] frac_dc_driver_pred;
    array[N_dc, n_sampling_times] simplex[2] frac_dc_clade_pred;
    array[N_clades_wt, n_sampling_times] simplex[2] frac_clade_wt_pred;


       for (t in 1:n_sampling_times) {

      // WT
      {
        real loc = fmin(1.0 - 1e-9, fmax(1e-9, frac_wt_lat[sampling_index[t]][1]));
        real p = beta_proportion_rng(loc, kappa);

        frac_wt_pred[t][1] = p;
        frac_wt_pred[t][2] = 1 - p;
      }


      // DRIVER
      for (i in 1:N_driver) {
        real loc = fmin(1.0 - 1e-9, fmax(1e-9, frac_driver_lat[i,sampling_index[t]][1]));
        real p = beta_proportion_rng(loc, kappa);

        frac_driver_pred[i,t][1] = p;
        frac_driver_pred[i,t][2] = 1 - p;
      }


      // DC DRIVER + CLADE
      for (i in 1:N_dc) {
        real loc1 = fmin(1.0 - 1e-9, fmax(1e-9, frac_driver_dc_lat[i,sampling_index[t]][1]));
        real p1 = beta_proportion_rng(loc1, kappa);

        frac_dc_driver_pred[i,t][1] = p1;
        frac_dc_driver_pred[i,t][2] = 1 - p1;

        real loc2 = fmin(1.0 - 1e-9, fmax(1e-9, frac_clade_dc_lat[i,t][1]));
        real p2 = beta_proportion_rng(loc2, kappa);

        frac_dc_clade_pred[i,t][1] = p2;
        frac_dc_clade_pred[i,t][2] = 1 - p2;
      }


      // CLADES WT
      for (i in 1:N_clades_wt) {
        real loc = fmin(1.0 - 1e-9, fmax(1e-9, frac_clade_wt_lat[i,t][1]));
        real p = beta_proportion_rng(loc, kappa);

        frac_clade_wt_pred[i,t][1] = p;
        frac_clade_wt_pred[i,t][2] = 1 - p;
      }
    }

    
    // =========================================================
    // COUNTS POSTERIOR PREDICTIVE
    // =========================================================
    
    array[N_driver] int m_driver_pred;
    array[N_driver_n] int m_driver_n_pred;
    array[N_dc,2] int m_dc_pred;
    array[N_clades_wt] int m_clade_wt_pred;
    int m_trunk_pred;
    
    m_trunk_pred = poisson_rng(2*mu*l*lambda_n*2*(tmrca-t_min));
    
      
        for (i in 1:N_driver) {
          
          m_driver_pred[i] =
          poisson_rng(
            4 * mu * l * lambda_n * (t_driver[i] - tmrca)
            );
        }
        
        for (i in 1:N_driver_n) {
          
          m_driver_n_pred[i] =
          poisson_rng(
            4 * mu * l * lambda_n * (t_driver_n[i] - tmrca)
            );
        }
        
        
        for (i in 1:N_dc) {
          
          m_dc_pred[i,1] =
          poisson_rng(
            2 * mu * l * lambda_n * 2 * (t_dc[i][1] - tmrca)
            );
            
            m_dc_pred[i,2] =
            poisson_rng(
              2 * mu * l * lambda_n * (1 + s_dc[i])
              * 2 * (t_dc[i][2] - t_dc[i][1])
              );
        }
        
        
        for (i in 1:N_clades_wt) {
          
          m_clade_wt_pred[i] =
          poisson_rng(
            2 * mu * l * lambda_n * 2 * (t_clade_wt[i] - tmrca)
            );
      }
      
          
  vector[N_driver_n] omega_p_driver_n_pred;
      
  for (i in 1:N_driver_n) {
    
    real W = exponential_rng(1);
    
    real U = uniform_rng(0,pi());
    
    real A = (sin((1.0 - alpha) * U) * pow(sin(alpha * U), alpha / (1.0 - alpha))) / pow(sin(U), 1.0 / (1.0 - alpha));
    
    real chi = pow(A / W, (1.0 - alpha) / alpha);
    
    real lambda_p = lambda_n * (1 + s_driver_n[i]) * (1 + s_epi);
    
    real zeta = 1 / (1 + s_epi);
    
    real target_ccf = 0.01;
    
    real sin_safe = sqrt(square(sin(pi() * zeta)) + 1e-8);
    
    real dt = times[sampling_index[n_sampling_times]] - t_driver_n[i];
    
    real Z2_target = target_ccf * zplus[n_sampling_times] / (chi + 1e-8);
    
    omega_p_driver_n_pred[i] = (lambda_p * sin_safe) / ( 
      pi() * bc_driver_n[i]) * pow( Z2_target * exp(-lambda_p * dt), zeta);
      
  }
  
   vector[N_driver_p] omega_n_driver_p_pred;
      
  for (i in 1:N_driver_p) {
     
     real lambda_p = lambda_n * (1 + s_driver_p[i]) * (1 + s_epi);
     
     real delta = lambda_p - lambda_n * (1 + s_driver_p[i]) + 1e-10;
     
     real dt = times[sampling_index[n_sampling_times]] - t_driver_p[i];
     
     real target_ccf = 0.01;
     
     real Z1_target = target_ccf * zminus[n_sampling_times];

     omega_n_driver_p_pred[i] = Z1_target*delta*exp(-lambda_p*dt);
        
  }
      
      
      
// ==============================
// POSTERIOR PREDICTIVE CCF (solo driver)
// ==============================

array[N_driver, n_sampling_times, 2] real ccf_driver_pred;
array[N_driver_n, n_sampling_times] real ccf_driver_n_pred;   
array[N_driver_p, n_sampling_times] real ccf_driver_p_pred;
array[N_dc, n_sampling_times, 2] real ccf_driver_dc_pred;
array[n_sampling_times] real frac_pos_pred;
array[n_sampling_times] real frac_neg_pred;


for (t in 1:n_sampling_times){

  real total_neg = Z_wt_lat[sampling_index[t]][1];
  real total_pos = Z_wt_lat[sampling_index[t]][2];

  for (i in 1:N_driver){
    total_neg += Z_driver_lat[i,sampling_index[t]][1];
    total_pos += Z_driver_lat[i,sampling_index[t]][2];
  }

  for (i in 1:N_driver_n){
    total_neg += exp(log_Z_driver_n_lat[i,sampling_index[t]]);
  }

for (i in 1:N_driver_p){
    total_pos += exp(log_Z_driver_p_lat[i,sampling_index[t]]);
  }
  
  for (i in 1:N_dc){
    total_neg += Z_driver_dc_lat[i,sampling_index[t]][1];
    total_pos += Z_driver_dc_lat[i,sampling_index[t]][2];
  }

  for (i in 1:N_driver){

    real p_neg = Z_driver_lat[i,sampling_index[t]][1] / total_neg;
    real p_pos = Z_driver_lat[i,sampling_index[t]][2] / total_pos;

    ccf_driver_pred[i,t][1] =
      beta_proportion_rng(p_neg, kappa);

    ccf_driver_pred[i,t][2] =
      beta_proportion_rng(p_pos, kappa);
  }

  for (i in 1:N_driver_n){

    real p = exp(log_Z_driver_n_lat[i,sampling_index[t]]) / total_neg;

    ccf_driver_n_pred[i,t] =
      beta_proportion_rng(p, kappa);
  }
  
    for (i in 1:N_driver_p){

    real p = exp(log_Z_driver_p_lat[i,sampling_index[t]]) / total_pos;

    ccf_driver_p_pred[i,t] =
      beta_proportion_rng(p, kappa);
  }

  for (i in 1:N_dc){

    real p_neg = Z_driver_dc_lat[i,sampling_index[t]][1] / total_neg;
    real p_pos = Z_driver_dc_lat[i,sampling_index[t]][2] / total_pos;

    ccf_driver_dc_pred[i,t][1] =
      beta_proportion_rng(p_neg, kappa);

    ccf_driver_dc_pred[i,t][2] =
      beta_proportion_rng(p_pos, kappa);
  }
  
frac_pos_pred[t] = beta_proportion_rng(total_pos / (total_neg + total_pos), kappa);
frac_neg_pred[t] = beta_proportion_rng(total_neg / (total_neg + total_pos), kappa);  
  
}

array[n_intermediate_times] real ztot_pred;

if (n_intermediate_times > 0) {

  for (k in 1:n_intermediate_times) {

    ztot_pred[k] = lognormal_rng(log(ztot_lat[k]), sigma_count);

  }

}
  
}


  
 

