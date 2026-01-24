

functions {
  vector switching_process(real t,
                          vector z,
                          array[] real theta //lambda_minus,lambda_plus, omega_minus, omega_plus
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

data {
  
  int n_times; 
  vector[5] z0;
  real t0;
  array[n_times] real zminus;
  array[n_times] real zplus;
  array[n_times] real t;
  real <lower = 0> alpha_lambda;
  real <lower = 0> beta_lambda;
  real ms_epi;
  real <lower=0> sigma_epi;
  real <lower = 0> alpha_n;
  real <lower = 0> beta_n;
  real <lower = 0> alpha_p;
  real <lower = 0> beta_p;
  // real <lower=0> alpha_minus;
  // real <lower=0> beta_minus;
  // real <lower=0> alpha_plus;
  // real <lower=0> beta_plus;
  
}

// transformed data {
//   real x_r[0];
//   int x_i[0];
// }

parameters {
  real<lower=0> lambda_minus;
  real<lower=0> s_hat;
  real<lower=0> omega_p; 
  real<lower=0> omega_n; 
  // real<lower=0> sigma_plus; 
  // real<lower=0> sigma_minus; 
}

transformed parameters {
  
  real <lower=-1> s_epi;
  s_epi = s_hat - 1;

  array[4] real theta;
  
  theta[1] = lambda_minus;
  theta[2] = lambda_minus*(1+s_epi);
  theta[3] = omega_n; // omega_minus
  theta[4] = omega_p; // omega_plus
  
  array[n_times] vector[5] z_hat = ode_rk45(switching_process, z0, t0, t, theta);
}

model {
  // real z_hat[n_times,5];
  
  target += gamma_lpdf(lambda_minus | alpha_lambda, beta_lambda);
  s_hat ~ lognormal(ms_epi,sigma_epi);
  target += gamma_lpdf(omega_p | alpha_p, beta_p);
  target += gamma_lpdf(omega_n | alpha_n, beta_n);
  
  // sigma_minus ~ gamma(alpha_minus,beta_minus);
  // sigma_plus ~ gamma(alpha_plus,beta_plus);

  // z_hat = integrate_ode_rk45(switching_process, z0, t0, t, theta, x_r, x_i);


  for (i in 1:n_times) {
    target += normal_lpdf(zminus[i] | z_hat[i,1], sqrt(z_hat[i,3]));
    target += normal_lpdf(zplus[i] | z_hat[i,2], sqrt(z_hat[i,4]));
  }

}

generated quantities {
  array [n_times] real pred_minus;
  array [n_times] real pred_plus;
  // real var_minus[n_times];
  // real var_plus[n_times];
  // real pred[n_times, 5] = integrate_ode_rk45(switching_process, z0, t0, t, theta, x_r, x_i);
  // array[n_times] vector[5] pred = ode_rk45(switching_process, z0, t0, t, theta);

  for (i in 1:n_times){
    pred_minus[i] = normal_rng(z_hat[i,1], sqrt(z_hat[i,3]));
    pred_plus[i] = normal_rng( z_hat[i,2], sqrt(z_hat[i,4]));
    // var_minus[i] = pred[i,3];
    // var_plus[i] = pred[i,4];
  }
  
  
  real  s_epi_prior;
  real <lower = 0> lambda_minus_prior;
  real <lower = 0> omega_p_prior;
  real <lower = 0> omega_n_prior;
  
  s_epi_prior = lognormal_rng(ms_epi,sigma_epi) - 1;
  lambda_minus_prior = gamma_rng(alpha_lambda,beta_lambda);
  omega_n_prior = gamma_rng(alpha_n,beta_n);
  omega_p_prior = gamma_rng(alpha_p,beta_p);
  
}


