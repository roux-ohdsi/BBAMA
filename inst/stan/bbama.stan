functions {
    // to convert a real to a (positive) integer
    int count_up(real how_far){
        int j = 0;
        while (j <= how_far)
        j = j + 1;
        return(j);
    }

    // interpolate from a sequence of values from log-likelihood function
    // doesn't deal with values outside the grid (proposal will be rejected with warning)
    // assumes they are ordered, with no duplicates
    real grid_lp(vector theta_tilde, vector[] ll, vector vals_evaled, int n_vals,
    real min_val, real max_val) {

        real lprob = 0; // collect the sum of the log-likelihoods

        for (i in 1:num_elements(theta_tilde)) {
            real x; // theta_tilde proposal
            real how_far; // index-ish of x in equally spaced vals_evaled vector
            int floor_ind; // index as integer
            int ceil_ind; // index of value above x
            real val_lo; // actual values evaluated at those indices
            real val_up;
            real ll_lo; // log-likelihoods at those indices
            real ll_up;

            x = theta_tilde[i];
            if (x < min_val || x > max_val)
              reject("proposal not in range");

            // find the index of the value right below x
            how_far = (x-min_val)*n_vals/(max_val-min_val);
            ceil_ind = count_up(how_far);
            if (ceil_ind == 1)
              reject("proposal not in range");
            floor_ind = ceil_ind - 1;
            // LL at loside grid value
            ll_lo = ll[floor_ind][i];
            // LL at upside grid value
            ll_up = ll[ceil_ind][i];
            val_lo = vals_evaled[floor_ind];
            val_up = vals_evaled[ceil_ind];
            // estimate LL at x based on this linear slope within interval around x
            lprob = ll_lo + (x - val_lo)*(ll_up - ll_lo)/(val_up - val_lo) + lprob;
        }

        return lprob;
    }
}

data {
    int<lower=1,upper=3> dist_main; // 1 = poisson, 2 = normal, 3 = grid
    int<lower=0,upper=3> dist_bias; // 0 = none, 1 = poisson, 2 = normal, 3 = grid

    int<lower=0,upper=1> fixed_main;
    int<lower=0,upper=1> fixed_bias;

    int<lower=1> M; // number of sites
    int<lower=0> N; // number of negative controls (overall)
    int<lower=0,upper=N> site[N]; // which site a negative control belongs to

    // poisson
    int y_int[dist_main == 1 ? M : 0]; // y_i for poisson
    int x_int[dist_bias == 1 ? N : 0]; // x_ij for poisson
    vector<lower=0>[dist_main == 1 ? M : 0] y_star; // y*
    vector<lower=0>[dist_bias == 1 ? N : 0] x_star; // x*_ij

    /// normal
    vector[dist_main == 2 ? M : 0] y; // log RRs
    vector[dist_bias == 2 ? N : 0] x; // log RRs of negative controls
    vector<lower=0>[dist_main == 2 ? M : 0] s_0; // standard errors of main effects
    vector<lower=0>[dist_bias == 2 ? N : 0] s_j; // standard errors of negative controls

    // grid
    int<lower=0> n_vals; // number of likelihood values per site (grid size)
    vector[dist_main == 3 ? M : 0] ll[n_vals]; // likelihood values (array of size M containing vectors with n_vals elements)
    vector[n_vals] vals_evaled; // values at which likelihood was evaluated
    real min_val; // min & max values at which likelihood evaluated
    real max_val;

    // fixed effects for main effect
    int K; // number of predictors
    matrix[fixed_main ? M : 0,K] p; // predictors

    // fixed effects for mean bias
    int L; // number of predictors
    matrix[fixed_bias ? M : 0,L] q; // predictors

    // prior distributions
    real mu_mean;
    real<lower=0> mu_sd;
    real tau_mean;
    real<lower=0> tau_sd;
    real lambda_mean;
    real<lower=0> lambda_sd;
    real eta_mean;
    real<lower=0> eta_sd;
    // can have separate priors for the sd
    // of the site-specific bias if you want
    vector[M] gamma_mean;
    vector<lower=0>[M] gamma_sd;
}

parameters {
    // parameters for the overall effect distribution
    real mu;
    real<lower=0> tau;

    // parameters for the overall bias distribution
    vector[dist_bias == 0 ? 0 : 1] lambda;
    vector<lower=0>[dist_bias == 0 ? 0 : 1] eta;

    //data source-specific parameters for the bias distribution
    vector[dist_bias == 0 ? 0 : M] delta;
    vector<lower=0>[dist_bias == 0 ? 0 : M] gamma;

    // true data source-specific effect
    vector[M] theta_0;
    // fixed effects for main effect
    vector[fixed_main ? K : 0] alpha;
    // fixed effects for mean bias
    vector[fixed_bias ? L : 0] omega;

    // true biases for the effects of interest and for negative controls
    vector[dist_bias == 0 ? 0 : M] beta_0;
    vector[N] beta_j;
}

model{
    vector[M] theta_tilde;
    vector[dist_bias == 0 ? 0 : M] delta_star;

    //poisson
    vector[dist_bias == 1 ? N : 0] beta_j_star;
    vector[dist_main == 1 ? M : 0] theta_tilde_star;

    // priors for overall effect
    mu ~ normal(mu_mean, mu_sd);
    tau ~ normal(tau_mean, tau_sd);

    // priors for bias distribution
    if (dist_bias > 0){
        lambda ~ normal(lambda_mean, lambda_sd);
        eta ~ normal(eta_mean, eta_sd);
        // data-source specific bias distribution
        for (i in 1:M) { // because I made lambda and eta vectors so they could be zero-length...
            delta[i] ~ normal(lambda, eta);
        }
        gamma ~ normal(gamma_mean, gamma_sd);

        if (fixed_bias == 1) {
            delta_star = delta + q * omega;
        } else {
            delta_star = delta;
        }

        //poisson
        if (dist_bias == 1) {
            for (j in 1:N){
                beta_j[j] ~ normal(delta_star[site[j]], gamma[site[j]]);
                beta_j_star[j] = exp(beta_j[j]) * x_star[j];
            }
            x_int ~ poisson(beta_j_star);
        } else {// normal
        for (j in 1:N){
            beta_j[j] ~ normal(delta_star[site[j]], gamma[site[j]]);
        }
        x ~ normal(beta_j, s_j);
        }
        beta_0 ~ normal(delta_star, gamma);
        // effect of interest
        if (fixed_main == 1) {
            theta_tilde = theta_0 + beta_0 + p * alpha;
        } else {
            theta_tilde = theta_0 + beta_0;
        }
    } else { // dist_bias == 0
        theta_0 ~ normal(mu, tau);
        if (fixed_main == 1) {
            theta_tilde = theta_0 + p * alpha;
        } else {
            theta_tilde = theta_0;
        }
    }

    //poisson
    if (dist_main == 1) {
        for (i in 1:M){
            theta_tilde_star[i] = exp(theta_tilde[i]) * y_star[i];
        }
        y_int ~ poisson(theta_tilde_star);
    } else if (dist_main == 2) {
        //normal
        y ~ normal(theta_tilde, s_0);
    } else {
        //grid
        target += grid_lp(theta_tilde, ll, vals_evaled, n_vals, min_val, max_val);
    }
}
generated quantities {
    vector[fixed_main ? M : 0] theta_transformed;
    vector[fixed_bias ? M : 0] delta_transformed;
    real<lower=0> tau_squared = tau^2;
    if (fixed_main) {
        theta_transformed = theta_0 + p * alpha;
    }
    if (fixed_bias) {
        delta_transformed = delta + q * omega;
    }
}
