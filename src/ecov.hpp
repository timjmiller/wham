template<class Type>
vector<Type> get_nll_Ecov(vector<int> Ecov_model, matrix<Type> Ecov_process_pars, matrix<Type> Ecov_re, vector<int> Ecov_use_re){
  /* 
     get nll contribtions for any Ecovs
       Ecov_model: which time series model to use (RW, AR1)
        Ecov_process_pars: parameters for time series models for each Ecov (columns)
          Ecov_re: random effects representing latent covariates
      Ecov_use_re: whether to include random effects in the likelihood
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_Ecov = Ecov_model.size();
  int n_y = Ecov_re.rows();
  vector<Type> nll_Ecov(n_Ecov); // nll contribution each Ecov_re
  nll_Ecov.setZero();

  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    // Ecov model option 1: RW
    if((Ecov_model(i) == 1) & (Ecov_use_re(i)==1)){
      Type Ecov_sig = exp(Ecov_process_pars(1,i)); // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
      Type Ecov1 = Ecov_process_pars(0,i);; // Ecov_x in year 1 (fixed effect)
      nll_Ecov(i) -= dnorm(Ecov_re(1,i), Ecov1, Ecov_sig, 1); // Ecov_re(0,i) set to NA
      for(int y = 2; y < n_y; y++){
        nll_Ecov(i) -= dnorm(Ecov_re(y,i), Ecov_re(y-1,i), Ecov_sig, 1);
      }
    }

    // Ecov model option 2: AR1
    if((Ecov_model(i) == 2) & (Ecov_use_re(i)==1)){
      //Type Ecov_mu = Ecov_process_pars(0,i);  // mean
      Type Ecov_phi = geninvlogit(Ecov_process_pars(2,i), Type(-1), Type(1), Type(1)); // autocorrelation
      //Type Ecov_phi = -1 + 2/(1 + exp(-Ecov_process_pars(2,i))); // autocorrelation
      Type Ecov_sig = exp(Ecov_process_pars(1,i)); // marginal sd
      vector<Type> re_i = Ecov_re.col(i);
      nll_Ecov(i) += SCALE(AR1(Ecov_phi), Ecov_sig)(re_i);
    }
  }
  return(nll_Ecov);
}

template<class Type>
matrix<Type> get_Ecov(vector<int> Ecov_model, matrix<Type> Ecov_process_pars, matrix<Type> Ecov_re, vector<int> Ecov_use_re){
  /* 
     get Ecov time series spanning years of environmental data and population model
       Ecov_model: which time series model to use (RW, AR1)
        Ecov_process_pars: parameters for time series models for each Ecov (columns)
          Ecov_re: random effects representing latent covariates
      Ecov_use_re: whether to include random effects in the likelihood
  */

  int n_Ecov = Ecov_model.size();
  int n_y = Ecov_re.rows();
  matrix<Type> Ecov_x(n_y, n_Ecov); // 'true' estimated Ecov (x_t in Miller et al. 2016 CJFAS)
  Ecov_x.setZero();

  // Ecov_model == 0) no Ecov
  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    // Ecov model option 1: RW
    if((Ecov_model(i) == 1) & (Ecov_use_re(i)==1)){
      Type Ecov1 = Ecov_process_pars(0,i); // Ecov_x in year 1 (fixed effect)
      Ecov_x(0,i) = Ecov1;
      for(int y = 1; y < n_y; y++) Ecov_x(1,i) = Ecov_re(1,i);
    }

    // Ecov model option 2: AR1
    if((Ecov_model(i) == 2) & (Ecov_use_re(i)==1)){
      Type Ecov_mu = Ecov_process_pars(0,i); // mean
      for(int y = 0; y < n_y; y++) Ecov_x(y,i) = Ecov_mu + Ecov_re(y,i);
    }
  } // end loop over Ecovs
  return(Ecov_x);
}

template<class Type>
matrix<Type> simulate_Ecov_re(vector<int> Ecov_model, matrix<Type> Ecov_process_pars, matrix<Type> Ecov_re, vector<int> Ecov_use_re){
  /* 
     simulate random effects for Ecov time series
       Ecov_model: which time series model to use (RW, AR1)
        Ecov_process_pars: parameters for time series models for each Ecov (columns)
          Ecov_re: random effects representing latent covariates
      Ecov_use_re: whether to simulate random effects (and whether to include likelihood contributions)
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_Ecov = Ecov_model.size();
  int n_y = Ecov_re.rows();
  matrix<Type> re_sim = Ecov_re;

  for(int i = 0; i < n_Ecov; i++){ // loop over Ecovs
    if(Ecov_use_re(i)==1){
      // Ecov model option 1: RW
      if(Ecov_model(i) == 1){
        Type Ecov1 = Ecov_process_pars(0,i);; // Ecov_x in year 1 (fixed effect)
        Type Ecov_sig = exp(Ecov_process_pars(1,i)); // sd (sig_x in Eq1, pg 1262, Miller et al. 2016)
        re_sim(1,i) = rnorm(Ecov1, Ecov_sig);
        for(int y = 2; y < n_y; y++){
          re_sim(y,i) = rnorm(re_sim(y-1,i), Ecov_sig);
        }
      }
      // Ecov model option 2: AR1
      if(Ecov_model(i) == 2){
        Type Ecov_phi = geninvlogit(Ecov_process_pars(2,i), Type(-1), Type(1), Type(1)); // autocorrelation
        //Type Ecov_phi = -1 + 2/(1 + exp(-Ecov_process_pars(2,i))); // autocorrelation
        Type Ecov_sig = exp(Ecov_process_pars(1,i)); // marginal sd
        vector<Type> re_i = re_sim.col(i);
        AR1(Ecov_phi).simulate(re_i);
        re_i *= Ecov_sig;
        for(int j = 0; j < re_i.size(); j++) re_sim(j,i) = re_i(j);
      }
    }
  }
  return(re_sim);
}

template<class Type>
matrix<Type> get_Ecov_out(matrix<Type> Ecov_x, int n_years_model, int n_years_proj, vector<int> ind_Ecov_out_start, vector<int> ind_Ecov_out_end){
  /* 
     get Ecov time series for years pertinent to specific effects on population
                    Ecov_x: Ecov time series spanning years of environmental data and population model 
             n_years_model: number of years for population model
              n_years_proj: number of years population model is projected
        ind_Ecov_out_start: which Ecov_x year is first year needed for effect on population
          ind_Ecov_out_end: which Ecov_x years is last year needed for effect on population
  */
  //int n_effects = Ecov_beta.dim(0); // 2 + n_indices (recruitment, mortality and any catchabilities)
  int n_Ecov = Ecov_x.cols();
  matrix<Type> Ecov_out(n_years_model + n_years_proj, n_Ecov); // Pop model uses Ecov_out(t) for processes in year t (Ecov_x shifted by lag and padded)
  Ecov_out.setZero(); // set Ecov_out = 0
  for(int i = 0; i < n_Ecov; i++){
      int ct = 0;
      for(int y = ind_Ecov_out_start(i); y < ind_Ecov_out_end(i) + 1 + n_years_proj; y++){
        Ecov_out(ct,i) = Ecov_x(y,i);
        ct++;
      }
  }
  return(Ecov_out);
}

template<class Type>
matrix<Type> get_Ecov_lm(matrix<Type>Ecov_beta, matrix<Type>Ecov_out, int n_years_model, int n_years_proj){
  /*
    Calculate ecov link model (possibly othogonal polynomial: b1*ecov + b2*ecov^2 + ...) --------------------
    to be added to linear model for particular parameter effects (R, q, M, mu, etc.).
    Ecov_beta: matrix (a sub-component of ecov_beta_*): n_ecov x n_poly
    Ecov_out: matrix returned by get_Ecov_out with Ecov time series spanning pertinent years
  */
  int n_Ecov = Ecov_out.cols();
  int n_poly = Ecov_beta.cols();
  // Ecov_lm stores the linear models for each Ecov each year.
  int ny = Ecov_out.rows(); //n_years_model + n_years_proj 
  matrix<Type> Ecov_lm(ny,n_Ecov); 
  Ecov_lm.setZero();
  for(int i = 0; i < n_Ecov; i++){
    vector<Type> thecol = Ecov_out.col(i);
    matrix<Type> X_poly(ny, n_poly);
    X_poly.setZero();
    if(n_poly == 1){ // n_poly = 1 if ecov effect is none or linear
      X_poly = thecol.matrix();
    } else { // n_poly > 1, get poly transformation for ith ecov
      X_poly = poly_trans(thecol, n_poly, n_years_model, n_years_proj); //need to treat X in projections differently
    }
    for(int y = 0; y < ny; y++){
      for(int j = 0; j < n_poly; j++){
        Ecov_lm(y,i) += Ecov_beta(i,j) * X_poly(y,j); // poly transformation returns design matrix, don't need to take powers
      }
    }
  }
  return Ecov_lm;
}

