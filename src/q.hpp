// Survey catchability functions

//prior on q
template<class Type>
vector<Type> get_nll_q_prior(vector<Type> q_prior_re, vector<Type> logit_q,  vector<Type> logit_q_prior_sigma, 
  vector<int> use_q_prior){
  /* 
    get any nll components for priors/posteriors on logit(catchability).
               q_prior_re: n_indices. RE for posterior (given prior) (mean) logit(catchability)
                  logit_q: n_indices. (mean) logit(catchability)
      logit_q_prior_sigma: n_indices. sd for prior distribution for RE logit(catchability)
              use_q_prior: n_indices. 0/1 whether to apply prior for each index
  */
  int n_indices = q_prior_re.size();
  vector<Type> nll_q_prior(n_indices);
  nll_q_prior.setZero();
  for(int i = 0; i < n_indices; i++) if(use_q_prior(i) == 1){
    //use prior for q? q_prior_re are random effects with mean logit_q (fixed) and sd = logit_q_prior_sigma.
    nll_q_prior(i) -= dnorm(q_prior_re(i), logit_q(i), logit_q_prior_sigma(i), 1);
  }
  return(nll_q_prior);
}
//done

//simulate from prior on q
template<class Type>
vector<Type> simulate_q_prior_re(vector<Type> q_prior_re, vector<Type> logit_q, 
  vector<Type> logit_q_prior_sigma, vector<int> use_q_prior){
  /* 
    simulate and RE for prior/posterior on logit(catchability).
               q_prior_re: n_indices. RE for posterior (given prior) (mean) logit(catchability)
                  logit_q: n_indices. (mean) logit(catchability)
      logit_q_prior_sigma: n_indices. sd for prior distribution for RE logit(catchability)
              use_q_prior: n_indices. 0/1 whether to apply prior for each index
  */
  
  int n_indices = q_prior_re.size();
  vector<Type> sim_q_prior_re(n_indices);
  sim_q_prior_re.setZero();
  for(int i = 0; i < n_indices; i++) if(use_q_prior(i) == 1){
    sim_q_prior_re(i) = rnorm(logit_q(i), logit_q_prior_sigma(i));
  }
  return(sim_q_prior_re);
}
//done

//auto-regressive random effects for q
template<class Type>
matrix<Type> get_nll_q(matrix<Type>q_repars, matrix<Type> q_re, vector<int> use_q_re){
  /* 
    get any nll components for time/age varying RE for logit(catchability).
      q_repars: n_indices x 2. parameters for distributions of random effects (sig, rho_y)
          q_re: n_years x n_indices. RE for logit(catchability).
      use_q_re: 0/1 whether to use time-varying random effects for each index
  */
  
  int n_indices = logit_q.size();
  int n_y = q_re.dim(0);
  matrix<Type> nll(n_y,n_indices);
  nll.setZero();
  vector<Type> sigma_q(n_indices);
  sigma_q.setZero();
  vector<Type> rho_q(n_indices);
  rho_q.setZero();
  for(int i = 0; i < n_indices; i++) {
    if(use_q_re(i) > 0) // random effects on q, q_re = AR1 deviations on (year,age), dim = n_years x n_M_a
    {
      sigma_q(i) = exp(q_repars(i,0)); // conditional sd
      rho_q(i) = rho_trans(q_repars(i,1)); // autocorrelation
      nll_q(0,i) -= dnorm(q_re(0,i), Type(0), sigma_q(i)*exp(-0.5 * log(1 - pow(rho_q(i),Type(2)))), 1);
      for(int y = 1; y < n_y; y++)
      {
        nll_q(y,i) -= dnorm(q_re(y,i), rho_q(i) * q_re(y-1,i), sigma_q(i), 1);
      }
    }
  }
  return(nll);
}
//done

//simulate auto-regressive random effects for q
template<class Type>
matrix<Type> simulate_q_re(matrix<Type>q_repars, matrix<Type> q_re, vector<int> use_q_re){
  /* 
    simulate andy time/age varying RE for logit(catchability).
      q_repars: n_indices x 2. parameters for distributions of random effects (sig, rho_y)
          q_re: n_years x n_indices. RE for logit(catchability).
      use_q_re: 0/1 whether to use time-varying random effects for each index
  */
  
  int n_indices = q_re.dim(1);
  int n_y = q_re.dim(0);
  matrix<Type> sim_q_re(n_y,n_indices);
  sim_q_re.setZero();
  vector<Type> sigma_q(n_indices);
  sigma_q.setZero();
  vector<Type> rho_q(n_indices);
  rho_q.setZero();
  for(int i = 0; i < n_indices; i++) {
    if(use_q_re(i) > 0) // random effects on q, q_re = AR1 deviations on year, dim = n_years x n_indices
    {
      sigma_q(i) = exp(q_repars(i,0)); // conditional sd
      rho_q(i) = rho_trans(q_repars(i,1)); // autocorrelation
      q_re(0,i) = rnorm(Type(0), sigma_q(i)*exp(-0.5 * log(1 - pow(rho_q(i),Type(2)))));
      for(int y = 1; y < n_y; y++) q_re(y,i) = rnorm(rho_q(i) * q_re(y-1,i), sigma_q(i));
    }
  }
  return(q_re);
}
//done

template<class Type>
matrix<Type> get_logit_q_mat(vector<Type> logit_q, matrix<Type> q_re, vector<Type> q_prior_re,  vector<int> use_q_prior, 
  vector<int> use_q_re, array<Type> Ecov_lm, matrix<int> use_Ecov) {
  /* 
    Construct n_years x n_indices matrices of logit(catchability)
    currently continues random processes in any projection years!
          logit_q: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
             q_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       q_prior_re: n_stocks x n_ages x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
      use_q_prior: n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
           Ecov_lm: (n_indices, n_years_pop, n_Ecov) linear predictor for any Ecov effects on trans_mu_base
          use_Ecov: n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 values indicating to use effects on migration for each stock for each region (less 1).
  */
  
  int n_indices = logit_q.size();
  int n_y = q_re.dim(0);
  matrix<Type> logit_q_mat(n_y, n_indices);
  logit_q_mat.setZero();
  for(int i = 0; i < n_indices; i++) {
    //use prior for q? q_prior_re are random effects with mean logit_q (fixed) and sd = logit_q_prior_sigma.
    if(use_q_prior(i) == 1){ 
      for(int y = 0; y < n_y; y++) logit_q_mat(y,i) += q_prior_re(i);
    }
    else for(int y = 0; y < n_y; y++) logit_q_mat(y,i) += logit_q(i);
    
    if(use_q_re(i) > 0) // random effects on q, q_re = AR1 deviations on (year,age), dim = n_years x n_M_a
    {
      for(int y = 0; y < n_y; y++)
      {
        logit_q_mat(y,i) += q_re(y,i); //add in q random effects.
      }
    }
    for(int y = 0; y < n_y; y++) {
      for(int j=0; j < n_Ecov; j++){
        if(use_Ecov(j,i) == 1){ // if ecov i affects q and which index
          logit_q_mat(y,i) += Ecov_lm(i,y,j);
        }
      }
    }
  }
  return(logit_q_mat);
}
//done

template<class Type>
matrix<Type> get_q(matrix<Type>logit_q_mat, vector<Type> q_lower, vector<Type> q_upper){
  /* 
    Construct n_years x n_indices matrix of catchability
      logit_q_mat: n_years_pop x n_indices. matrix retruned by get_logit_q_mat
          q_lower: n_indices. Lower bounds of catchability (default = 0)
          q_upper: n_indices. Upper bounds of catchability (default = 1000)
  */
  int n_y = logit_q_mat.dim(0);
  int n_ind = logit_q_mat.dim(1);
  matrix<Type> q(n_y, n_ind);
  q.setZero();
  for(int y = 0; y < n_y; y++) {
    for(int ind = 0; ind < n_ind; ind++) {
      q(y,ind) = inv_logit(logit_q_mat(y,ind),q_lower(ind),q_upper(ind),1);
      //q(y,ind) = q_lower(ind) + (q_upper(ind) - q_lower(ind))/(1 + exp(-logit_q_mat(y,ind)));
    }
  }
  return(q);
}
//done

template<class Type>
array<Type> get_QAA(matrix<Type> q, vector<matrix<Type>> selAA, matrix<int> selblock_pointer, int n_years_model, int n_ages){
  /* 
    Construct n_indices x n_years x n_ages array of catchability at age (q*sel)
    NOTE: q may be projected but currently selectivity in last model year is used.
      q: n_years_pop x n_indices. matrix retruned by get_logit_q_mat
          selAA: n_selblocks vector of n_years x n_ages selectivity matrices
          selblock_pointer: n_years x n_indices. which selectivity at age to use from selAA for each year
          n_years_model
  */
  int n_y = q.dim(0);
  int n_ind = q.dim(1);
  array<Type> QAA(n_ind,n_y,n_ages)
  QAA.setZero();
  // Construct survey catchability-at-age (QAA)
  for(int i = 0; i < n_ind; i++)
  {
  // add ecov effect on M (by year, shared across ages)
    for(int y = 0; y < n_y; y++)
    {
      if(y < n_years_model) {
        for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(y,i) * selAA(selblock_pointer_indices(y,i)-1)(y,a);
      } else { //for projections, just use last years selectivity for now
        for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(y,i) * selAA(selblock_pointer_indices(n_years_model-1,i)-1)(n_years_model-1,a);
      }
    }
  }
  return(QAA);
}
//done
