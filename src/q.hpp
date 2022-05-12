// Survey catchability functions

//prior on q
template<class Type>
vector<Type> get_nll_q_prior(vector<Type> q_prior_re, vector<Type> logit_q, 
  vector<Type> logit_q_prior_sigma, vector<int> use_q_prior){
  int n_indices = q_prior_re.size();
  vector<Type> nll_q_prior(n_indices);
  nll_q_prior.setZero();
  for(int i = 0; i < n_indices; i++) if(use_q_prior(i) == 1){
    //use prior for q? q_prior_re are random effects with mean logit_q (fixed) and sd = logit_q_prior_sigma.
    nll_q_prior(i) -= dnorm(q_prior_re(i), logit_q(i), logit_q_prior_sigma(i), 1);
  }
  return(nll_q_prior);
}

//simulate from prior on q
template<class Type>
vector<Type> simulate_q_prior_re(vector<Type> q_prior_re, vector<Type> logit_q, 
  vector<Type> logit_q_prior_sigma, vector<int> use_q_prior){
  
  int n_indices = q_prior_re.size();
  vector<Type> sim_q_prior_re(n_indices);
  sim_q_prior_re.setZero();
  for(int i = 0; i < n_indices; i++) if(use_q_prior(i) == 1){
    sim_q_prior_re(i) = rnorm(logit_q(i), logit_q_prior_sigma(i));
  }
  return(sim_q_prior_re);
}

//auto-regressive random effects for q
template<class Type>
matrix<Type> get_nll_q(matrix<Type>q_repars, matrix<Type> q_re, vector<int> use_q_re){
  
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

//simulate auto-regressive random effects for q
template<class Type>
matrix<Type> simulate_q_re(matrix<Type>q_repars, matrix<Type> q_re, vector<int> use_q_re){
  
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

template<class Type>
matrix<Type> get_logit_q_mat(vector<Type> logit_q, matrix<Type> q_re, vector<Type> q_prior_re, array<Type> Ecov_lm_q, 
  vector<int> use_q_prior_re, vector<int> use_q_re){
  
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
        if(use_Ecov_q(j,i) == 1){ // if ecov i affects q and which index
          logit_q_mat(y,i) += Ecov_lm_q(j,i,y);
        }
      }
    }
  }
  return(logit_q_mat);
}

template<class Type>
matrix<Type> get_q(matrix<Type>logit_q_mat, vector<Type> q_lower, vector<Type> q_upper)
{
  int n_y = logit_q_mat.dim(0);
  int n_ind = logit_q_mat.dim(1);
  matrix<Type> q(n_y, n_ind);
  q.setZero();
  for(int y = 0; y < n_y; y++) {
    for(int ind = 0; ind < n_ind; ind++) {
      q(y,ind) = q_lower(ind) + (q_upper(ind) - q_lower(ind))/(1 + exp(-logit_q_mat(y,ind)));
    }
  }
  return(q);
}

template<class Type>
array<Type> get_QAA(matrix<Type> q, vector<matrix<Type>> selAA, matrix<int> selblock_pointer, vector<Type> q_lower, int n_years_model, int n_ages)
{  
  int n_y = q.dim(0);
  int n_ind = q.dim(1);
  array<Type> QAA(n_ind,n_y,n_ages)
  // Construct survey catchability-at-age (QAA)
  for(int i = 0; i < n_ind; i++)
  {
  // add ecov effect on M (by year, shared across ages)
    for(int y = 0; y < n_years_model; y++)
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(y,i) * selAA(selblock_pointer_indices(y,i)-1)(y,a);
    }
    //for projections, just use last years selectivity for now
    if(ny > n_years_model) for(int y = n_years_model; y < n_y; y++) for(int a = 0; a < n_ages; a++) {
      QAA(i,y,a) = q(y,i) * selAA(selblock_pointer_indices(n_years_model-1,i)-1)(n_years_model-1,a);
    }
  }
  return(QAA);
}
