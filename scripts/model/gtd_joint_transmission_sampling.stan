functions {
  real dirichlet_multinomial_custom_lpmf(int[] y, vector alpha) {
    real alpha_sum = sum(alpha);
    int n = sum(y);
    real lp = lgamma(n + 1) + lgamma(alpha_sum) - lgamma(n + alpha_sum);
    for (l in 1:num_elements(y))
      lp += lgamma(y[l] + alpha[l]) - lgamma(alpha[l]) - lgamma(y[l] + 1);
    return lp;
  }
}

data {
  int<lower=1> C;
  int<lower=2> T;
  int<lower=2> L;
  int<lower=1> P;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=0, upper=1> use_project_effects;
  int<lower=0> cases[C, T];
  matrix[T, K] B;
  matrix[C, T] reporting_change;
  matrix<lower=0>[C, L] initial_alpha;
  real<lower=0> import_exposure[C, T, L];
  int<lower=1, upper=C> obs_country[J];
  int<lower=1, upper=T> obs_month[J];
  int<lower=1, upper=P> obs_project[J];
  int<lower=0> y_genome[J, L];
}

parameters {
  matrix[C, K] r_coef;
  vector[L] log_theta_raw;
  vector<lower=0>[C] import_scale;
  vector<lower=0>[C] density_feedback;
  simplex[L] q0[C];
  matrix[P, L] project_raw;
  real<lower=0> sigma_project;
  vector[C] reporting_jump;
  vector<lower=0>[C] phi_cases;
  real<lower=0> genome_concentration;
}

transformed parameters {
  vector[L] log_theta;
  real<lower=0, upper=1> q[C, T, L];
  real<lower=0> mu_cases[C, T];
  real<lower=0> local_component[C, T, L];
  real<lower=0> import_component[C, T, L];
  matrix[P, L] project_effect;

  log_theta = log_theta_raw - mean(log_theta_raw);
  for (p in 1:P) {
    real row_mean = mean(to_vector(project_raw[p]));
    for (l in 1:L)
      project_effect[p, l] = use_project_effects * sigma_project *
                             (project_raw[p, l] - row_mean);
  }

  for (c in 1:C) {
    for (l in 1:L) {
      q[c, 1, l] = q0[c, l];
      local_component[c, 1, l] = 0;
      import_component[c, 1, l] = 0;
    }
    mu_cases[c, 1] = cases[c, 1] + 1e-6;

    for (t in 2:T) {
      real total_component = 0;
      real log_r = dot_product(to_vector(B[t]), to_vector(r_coef[c])) -
                   density_feedback[c] * log1p(cases[c, t - 1] / 1000.0);
      for (l in 1:L) {
        local_component[c, t, l] =
          exp(log_r + log_theta[l]) * q[c, t - 1, l] * (cases[c, t - 1] + 0.5);
        import_component[c, t, l] =
          import_scale[c] * import_exposure[c, t, l];
        total_component += local_component[c, t, l] +
                           import_component[c, t, l] + 1e-9;
      }
      for (l in 1:L)
        q[c, t, l] = (
          local_component[c, t, l] + import_component[c, t, l] + 1e-9
        ) / total_component;
      mu_cases[c, t] = total_component * exp(
        reporting_jump[c] *
        (reporting_change[c, t] - reporting_change[c, t - 1])
      );
    }
  }
}

model {
  to_vector(r_coef) ~ normal(0, 0.8);
  for (c in 1:C)
    for (k in 2:K)
      r_coef[c, k] - r_coef[c, k - 1] ~ normal(0, 0.35);

  log_theta_raw ~ normal(0, 0.20);
  import_scale ~ lognormal(log(100), 1.5);
  density_feedback ~ normal(0.5, 0.35);
  to_vector(project_raw) ~ normal(0, 1);
  sigma_project ~ normal(0, 0.5);
  reporting_jump ~ normal(0, 0.5);
  phi_cases ~ lognormal(log(10), 1);
  genome_concentration ~ gamma(2, 0.1);

  for (c in 1:C)
    q0[c] ~ dirichlet(to_vector(initial_alpha[c]));

  for (c in 1:C)
    for (t in 2:T)
      cases[c, t] ~ neg_binomial_2(mu_cases[c, t], phi_cases[c]);

  for (j in 1:J) {
    vector[L] log_p;
    vector[L] alpha;
    for (l in 1:L)
      log_p[l] = log(q[obs_country[j], obs_month[j], l]) +
                 project_effect[obs_project[j], l];
    alpha = genome_concentration * softmax(log_p);
    target += dirichlet_multinomial_custom_lpmf(y_genome[j] | alpha);
  }
}

generated quantities {
  vector[L] lineage_relative_transmission;
  real<lower=0, upper=1> post_import_fraction[C];
  int<lower=0> cases_rep[C, T];

  lineage_relative_transmission = exp(log_theta);

  for (c in 1:C) {
    real post_import = 0;
    real post_total = 0;
    cases_rep[c, 1] = cases[c, 1];
    for (t in 2:T) {
      if (t >= 49) {
        for (l in 1:L) {
          post_import += import_component[c, t, l];
          post_total += local_component[c, t, l] + import_component[c, t, l];
        }
      }
      cases_rep[c, t] = neg_binomial_2_rng(mu_cases[c, t], phi_cases[c]);
    }
    post_import_fraction[c] = post_total > 0 ? post_import / post_total : 0;
  }
}
