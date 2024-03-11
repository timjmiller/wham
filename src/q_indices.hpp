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
matrix<Type> get_nll_q_re(matrix<Type>q_repars, matrix<Type> q_re, vector<int> use_q_re, vector<int> years_use){
  /* 
    get any nll components for time/age varying RE for logit(catchability).
      q_repars: n_indices x 2. parameters for distributions of random effects (sig, rho_y)
          q_re: n_years x n_indices. RE for logit(catchability).
      use_q_re: 0/1 whether to use time-varying random effects for each index
     years_use: is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  
  int n_indices = q_re.cols();
  int n_y = years_use.size(); 
  // int n_y = q_re.rows();
  matrix<Type> nll_q(n_y,n_indices);
  nll_q.setZero();
  vector<Type> sigma_q(n_indices);
  sigma_q.setZero();
  vector<Type> rho_q(n_indices);
  rho_q.setZero();
  for(int i = 0; i < n_indices; i++) {
    if(use_q_re(i) > 0) // random effects on q, q_re = AR1 deviations on (year,age), dim = n_years x n_M_a
    {
      sigma_q(i) = exp(q_repars(i,0)); // conditional sd
      rho_q(i) = geninvlogit(q_repars(i,1),Type(-1), Type(1),Type(1)); // autocorrelation. using scale =1 ,2 is legacy
      nll_q(years_use(0),i) -= dnorm(q_re(years_use(0),i), Type(0), sigma_q(i)*exp(-0.5 * log(1 - pow(rho_q(i),Type(2)))), 1);
      for(int y = 1; y < n_y; y++)
      {
        nll_q(years_use(y),i) -= dnorm(q_re(years_use(y),i), rho_q(i) * q_re(years_use(y)-1,i), sigma_q(i), 1);
      }
    }
  }
  return(nll_q);
}
//done

//simulate auto-regressive random effects for q
template<class Type>
matrix<Type> simulate_q_re(matrix<Type>q_repars, matrix<Type> q_re, vector<int> use_q_re, vector<int> years_use){
  /* 
    simulate andy time/age varying RE for logit(catchability).
      q_repars: n_indices x 2. parameters for distributions of random effects (sig, rho_y)
          q_re: n_years x n_indices. RE for logit(catchability).
      use_q_re: 0/1 whether to use time-varying random effects for each index
     years_use: is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  
  int n_indices = q_re.cols();
  int n_y = years_use.size();
  matrix<Type> sim_q_re = q_re;
  //sim_q_re.setZero();
  vector<Type> sigma_q(n_indices);
  sigma_q.setZero();
  vector<Type> rho_q(n_indices);
  rho_q.setZero();
  for(int i = 0; i < n_indices; i++) {
    if(use_q_re(i) > 0) // random effects on q, q_re = AR1 deviations on year, dim = n_years x n_indices
    {
      sigma_q(i) = exp(q_repars(i,0)); // conditional sd
      rho_q(i) = geninvlogit(q_repars(i,1),Type(-1), Type(1),Type(1)); // autocorrelation, using scale =1 ,2 is legacy
      sim_q_re(years_use(0),i) = rnorm(Type(0), sigma_q(i)*exp(-0.5 * log(1 - pow(rho_q(i),Type(2)))));
      for(int y = 1; y < n_y; y++) sim_q_re(years_use(y),i) = rnorm(rho_q(i) * sim_q_re(years_use(y)-1,i), sigma_q(i));
    }
  }
  return(sim_q_re);
}
//done

template<class Type>
matrix<Type> get_logit_q_mat(vector<Type> logit_q, matrix<Type> q_re, vector<Type> q_prior_re,  vector<int> use_q_prior, 
  vector<int> use_q_re, array<Type> Ecov_lm, matrix<int> Ecov_how) {
  /* 
    Construct n_years x n_indices matrices of logit(catchability)
    currently continues random processes in any projection years!
          logit_q: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
             q_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       q_prior_re: n_stocks x n_ages x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
      use_q_prior: n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
           Ecov_lm: (n_indices, n_years_pop, n_Ecov) linear predictor for any Ecov effects on trans_mu_base
          Ecov_how: n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 values indicating to use effects on migration for each stock for each region (less 1).
  */
  
  int n_indices = q_re.cols();
  int n_y = q_re.rows();
  int n_Ecov = Ecov_lm.dim(2);
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
        if(Ecov_how(j,i) == 1){ // if ecov i affects q and which index
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
  int n_y = logit_q_mat.rows();
  int n_ind = logit_q_mat.cols();
  matrix<Type> q(n_y, n_ind);
  q.setZero();
  for(int y = 0; y < n_y; y++) {
    for(int ind = 0; ind < n_ind; ind++) {
      q(y,ind) = geninvlogit(logit_q_mat(y,ind),q_lower(ind),q_upper(ind),Type(1));
      //q(y,ind) = q_lower(ind) + (q_upper(ind) - q_lower(ind))/(1 + exp(-logit_q_mat(y,ind)));
    }
  }
  return(q);
}
//done

template<class Type>
array<Type> get_QAA(matrix<Type> q, vector<matrix<Type>> selAA, matrix<int> selblock_pointer_indices, int n_years_model, int n_ages){
  /* 
    Construct n_indices x n_years x n_ages array of catchability at age (q*sel)
    NOTE: q may be projected but currently selectivity in last model year is used.
      q: n_years_pop x n_indices. matrix retruned by get_logit_q_mat
          selAA: n_selblocks vector of n_years x n_ages selectivity matrices
          selblock_pointer: n_years x n_indices. which selectivity at age to use from selAA for each year
          n_years_model
  */
  int n_y = q.rows();
  int n_ind = q.cols();
  array<Type> QAA(n_ind,n_y,n_ages);
  QAA.setZero();
  // Construct survey catchability-at-age (QAA)
  for(int i = 0; i < n_ind; i++)
  {
  // add ecov effect on M (by year, shared across ages)
    for(int y = 0; y < n_y; y++)
    {
      if(y < n_years_model) {
        for(int a = 0; a < n_ages; a++) QAA(i,y,a) = q(y,i) * selAA(selblock_pointer_indices(y,i)-1)(y,a);
      } else { //for projections, just use last years selectivity for now
        for(int a = 0; a < n_ages; a++) QAA(i,y,a) = q(y,i) * selAA(selblock_pointer_indices(n_years_model-1,i)-1)(n_years_model-1,a);
      }
    }
  }
  return(QAA);
}
//done

template<class Type>
array<Type> get_pred_IAA(array<Type>QAA, array<Type> NAA_index) {

  int n_stocks = NAA_index.dim(0);
  int n_indices = NAA_index.dim(1);
  int n_y = NAA_index.dim(2); //n_years_model
  int n_ages = NAA_index.dim(3);
  array<Type> pred_IAA(n_indices, n_y,n_ages);
  pred_IAA.setZero();
  for(int s = 0; s < n_stocks; s++) for(int i = 0; i < n_indices; i++) for(int y = 0; y < n_y; y++) {
    for(int a = 0; a < n_ages; a++) {
      //get numbers at age a for stock s in each region at time of spawning
      pred_IAA(i,y,a) += QAA(i,y,a) * NAA_index(s,i,y,a);
    }
  }
  return pred_IAA;
}

template<class Type>
matrix<Type> get_pred_indices(array<Type> pred_IAA, vector<int> units_indices, array<Type> waa, vector<int> waa_pointer_indices){
  int n_indices = pred_IAA.dim(0);
  int n_y = pred_IAA.dim(1);
  int n_ages = pred_IAA.dim(2);
  matrix<Type> pred_indices(n_y, n_indices);
  pred_indices.setZero();

  for(int i = 0; i < n_indices; i++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) {
    if(units_indices(i) == 1) {
      pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(i,y,a);
    } else pred_indices(y,i) += pred_IAA(i,y,a);    
  }
  return pred_indices;
}

template <class Type>
matrix<Type> get_pred_log_indices(matrix<Type> pred_indices, matrix<Type> agg_index_sigma, vector<Type> log_index_sig_scale, 
  int bias_correct_oe){
  int n_y = agg_index_sigma.rows();
  int n_indices = agg_index_sigma.cols();
  matrix<Type> pred_log_indices(n_y,n_indices);
  pred_log_indices.setZero();

  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_indices; i++) {
    pred_log_indices(y,i) = log(pred_indices(y,i));
    Type sig = agg_index_sigma(y,i)*exp(log_index_sig_scale(i));
    if(bias_correct_oe) pred_log_indices(y,i) -= 0.5*exp(2*log(sig));
  }
  return pred_log_indices;
}

template <class Type>
matrix<Type> get_nll_agg_indices(matrix<Type> pred_log_indices, matrix<Type> agg_index_sigma, vector<Type> log_index_sig_scale,
  vector<Type> obsvec, matrix<int> use_indices, matrix<int> keep_I, data_indicator<vector<Type>, Type> keep){
  int n_y = agg_index_sigma.rows();
  int n_indices = agg_index_sigma.cols();

  matrix<Type> nll_agg_indices(n_y,n_indices);
  nll_agg_indices.setZero();

  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_indices; i++) if(use_indices(y,i)){
    Type sig = agg_index_sigma(y,i)*exp(log_index_sig_scale(i));
    nll_agg_indices(y,i) -= keep(keep_I(y,i)) * dnorm(obsvec(keep_I(y,i)), pred_log_indices(y,i), sig,1);
    nll_agg_indices(y,i) -= keep.cdf_lower(keep_I(y,i)) * log(squeeze(pnorm(obsvec(keep_I(y,i)), pred_log_indices(y,i), sig)));
    nll_agg_indices(y,i) -= keep.cdf_upper(keep_I(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_I(y,i)), pred_log_indices(y,i), sig)));
  }
  return nll_agg_indices;
}

template <class Type>
matrix<Type> simulate_agg_indices(matrix<Type> pred_log_indices, matrix<Type> agg_indices, matrix<Type> agg_index_sigma, 
  vector<Type> log_index_sig_scale, matrix<int> use_indices){
  int n_y = agg_indices.rows();
  int n_indices = agg_indices.cols();
  matrix<Type> agg_indices_out = agg_indices;
  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_indices; i++) if(use_indices(y,i)) {
    Type sig = agg_index_sigma(y,i)*exp(log_index_sig_scale(i));
    agg_indices(y,i) = exp(rnorm(pred_log_indices(y,i), sig));
  }
  return agg_indices;
}

template <class Type>
vector<Type> sim_agg_indices_in_obsvec(vector<Type> obsvec, matrix<int> keep_I, matrix<Type> agg_indices, matrix<int> use_indices){
  int n_y = use_indices.rows();
  int n_indices = use_indices.cols();
  vector<Type> obsvec_out = obsvec;
  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_indices; i++){
    if(use_indices(y,i)) obsvec_out(keep_I(y,i)) = log(agg_indices(y,i));
  }
  return obsvec_out;
}

template <class Type>
array<Type> get_pred_index_paa(array<Type> pred_IAA, vector<int> units_index_paa, 
  array<Type> waa, vector<int> waa_pointer_indices){
  int n_indices = pred_IAA.dim(0);
  int n_y = pred_IAA.dim(1);
  int n_ages = pred_IAA.dim(2);
  array<Type> pred_index_paa(n_indices,n_y, n_ages);

  for(int i = 0; i < n_indices; i++) for(int y = 0; y < n_y; y++){
    Type tsum = 0.0;
    for(int a = 0; a < n_ages; a++){
      if(units_index_paa(i) == 1) pred_IAA(i,y,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(i,y,a);
      tsum += pred_IAA(i,y,a);
    }
    for(int a = 0; a < n_ages; a++){
      pred_index_paa(i,y,a) = pred_IAA(i,y,a)/tsum;
    }
  }
  return pred_index_paa;
}

template <class Type>
matrix<Type> get_nll_index_acomp(array<Type> pred_index_paa, matrix<int> use_index_paa, array<Type> index_paa,
  matrix<Type> index_Neff, vector<int> age_comp_model_indices, matrix<Type> index_paa_pars, 
  array<int> keep_Ipaa, data_indicator<vector<Type>, Type> keep, vector<Type> obsvec, vector<int> agesvec, int do_osa){
  int n_indices = pred_index_paa.dim(0);
  int n_y = index_paa.dim(1);
  int n_ages = pred_index_paa.dim(2);
  matrix<Type> nll_index_acomp(n_y,n_indices);
  nll_index_acomp.setZero();

  for(int i = 0; i < n_indices; i++) for(int y = 0; y < n_y; y++)if(use_index_paa(y,i)) {
    vector<Type> paa_obs_y(n_ages);
    vector<Type> t_pred_paa(n_ages);
    for(int a = 0; a < n_ages; a++){
      t_pred_paa(a) = pred_index_paa(i,y,a);
      paa_obs_y(a) = index_paa(i,y,a);
    }
    //NB: indexing in obsvec MUST be: keep_Ipaa(i,y,0),...,keep_Ipaa(i,y,0) + keep_Ipaa(i,y,1) - 1
    //keep_Ipaa(i,y,0) is first val, keep_Ipaa(i,y,1) is the length of the vector
    vector<Type> tf_paa_obs = obsvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
    vector<int> ages_obs_y = agesvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
    nll_index_acomp(y,i) -= get_acomp_ll(tf_paa_obs, t_pred_paa, index_Neff(y,i), ages_obs_y, age_comp_model_indices(i), 
      vector<Type>(index_paa_pars.row(i)), keep.segment(keep_Ipaa(i,y,0),keep_Ipaa(i,y,1)), do_osa, paa_obs_y);
  }
  return nll_index_acomp;
}

template <class Type>
vector<Type> simulate_index_paa_in_obsvec(vector<Type> obsvec, vector<int> agesvec, array<Type> pred_index_paa, matrix<int> use_index_paa,
  array<int> keep_Ipaa, matrix<Type> index_Neff, vector<int> age_comp_model_indices, matrix<Type> index_paa_pars){
  int n_indices = pred_index_paa.dim(0);
  int n_y = keep_Ipaa.dim(1);
  int n_ages = pred_index_paa.dim(2);
  vector<Type> obsvec_out = obsvec;
  for(int i = 0; i < n_indices; i++) for(int y = 0; y < n_y; y++) if(use_index_paa(y,i)) {
    vector<Type> t_pred_paa(n_ages);
    for(int a = 0; a < n_ages; a++) t_pred_paa(a) = pred_index_paa(i,y,a);
    vector<int> ages_obs_y = agesvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
    vector<Type> tf_paa_obs = sim_acomp(t_pred_paa, index_Neff(y,i), ages_obs_y, age_comp_model_indices(i), 
      vector<Type>(index_paa_pars.row(i)));
    obsvec_out.segment(keep_Ipaa(i,y,0),keep_Ipaa(i,y,1)) = tf_paa_obs;
  }
  return obsvec_out;
}

template <class Type>
array<Type> sim_obsvec_in_index_paa(vector<Type> obsvec, vector<int> agesvec, array<Type> index_paa, matrix<int> use_index_paa, array<int> keep_Ipaa, 
  vector<int> age_comp_model_indices){
  int n_indices = index_paa.dim(0);
  int n_y = index_paa.dim(1);
  int n_ages = index_paa.dim(2);
  array<Type> index_paa_out = index_paa;
  vector<Type> paa_obs_y(n_ages);
  for(int i = 0; i < n_indices; i++) for(int y = 0; y < n_y; y++) if(use_index_paa(y,i)) {
    vector<int> ages_obs_y = agesvec.segment(keep_Ipaa(i,y,0), keep_Ipaa(i,y,1));
    vector<Type> tf_paa_obs = obsvec.segment(keep_Ipaa(i,y,0),keep_Ipaa(i,y,1));
    paa_obs_y = make_paa(tf_paa_obs, age_comp_model_indices(i), ages_obs_y, n_ages);
    for(int a = 0; a < n_ages; a++) index_paa_out(i,y,a) = paa_obs_y(a);
  }
  return index_paa_out;
}

template <class Type>
vector<Type> obsvec_to_paa(int i, int y, vector<Type> obsvec, vector<int> agesvec, matrix<int> use_paa, array<int> keep_paa, 
  vector<int> age_comp_model, int n_ages){
  vector<Type> paa(n_ages);
  paa.setZero();
  if(use_paa(y,i)) {
    vector<int> ages_obs_y = agesvec.segment(keep_paa(i,y,0), keep_paa(i,y,1));
    vector<Type> tf_paa_obs = obsvec.segment(keep_paa(i,y,0),keep_paa(i,y,1));
    paa = make_paa(tf_paa_obs, age_comp_model(i), ages_obs_y, paa);
  }
  return paa;
}
