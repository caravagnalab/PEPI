

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


  for (i in 1:n_times) {
    target += normal_lpdf(zminus[i] | z_hat[i,1], sqrt(z_hat[i,3]));
    target += normal_lpdf(zplus[i] | z_hat[i,2], sqrt(z_hat[i,4]));
  }

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


