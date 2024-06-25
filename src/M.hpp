template<class Type>
matrix<Type> get_nll_M(array<Type> M_repars, matrix<int> M_re_model, int M_model, array<Type>M_re, matrix<int> n_M_re, vector<int> years_use){
  /* 
     get nll contribtions for any M random effects
          M_repars: sd and correlation parameters for M random effects
        M_re_model: (n_stocks x n_regions) 1: no RE, 2: iid/ar1 RE for age (constant by year), 3: iid/ar1 RE year (constant by age), 4: 2D iid/AR1 (age, year) (not yet 5: 3DAR1 (age,year,region))
           M_model: 1: M= f(age), 2: M = f(WAA)
              M_re: random effects for 
            n_M_re: (n_stocks x n_regions) number of random effects (e.g., number of age classes or number of M=f(WAA) pars) each year), max(n_M_re)<=n_ages)
         years_use: is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = M_re.dim(0);
  int n_regions = M_re.dim(1);
  int n_y = years_use.size();
  //int n_y = M_re.dim(2);
  matrix<Type> nll_M(n_stocks,n_regions);
  nll_M.setZero();
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r< n_regions; r++) {
    if(M_re_model(s,r) > 1){ // random effects on M, M_re = 2D AR1 deviations on M(year,age), dim = n_years x n_M_re(s,r)
      Type sigma_M = exp(M_repars(s,r,0));
      Type rho_M_a = geninvlogit(M_repars(s,r,1),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy
      Type rho_M_y = geninvlogit(M_repars(s,r,2),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy
      Type Sigma_M;
      // likelihood of M deviations, M_re
      array<Type> M_re_r_s(n_y,n_M_re(s,r));
      M_re_r_s.setZero();
      for(int y = 0; y < n_y; y++)for(int a = 0; a < n_M_re(s,r); a++) M_re_r_s(y,a) = M_re(s,r,years_use(y),a); //first n_M_re(s,r) columns
      if(M_re_model(s,r) == 4){ //2D AR1: age, year (could be iid)
        Sigma_M = pow(pow(sigma_M,2) / ((1-pow(rho_M_y,2))*(1-pow(rho_M_a,2))),0.5);
        nll_M(s,r) += SCALE(SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)), Sigma_M)(M_re_r_s); // must be array, not matrix!
      } else {
        if(M_re_model(s,r) == 2){ // 1D ar1_a
          vector<Type> Mre0 = M_re_r_s.matrix().row(0);
          Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_a,2)),0.5);
          nll_M(s,r) += SCALE(AR1(rho_M_a), Sigma_M)(Mre0);
        } else { // M_re_model(s,r) = 3, 1D ar1_y
          vector<Type> Mre0 = M_re_r_s.matrix().col(0); //just first column
          Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_y,2)),0.5);
          nll_M(s,r) += SCALE(AR1(rho_M_y), Sigma_M)(Mre0);
        }
      }
    }
  }
  return(nll_M);
}
//done

template<class Type>
array<Type> simulate_M_re(array<Type> M_repars, matrix<int> M_re_model, int M_model, array<Type>M_re, matrix<int> n_M_re, vector<int> years_use){
  /* 
     simulate M random effects
          M_repars: sd and correlation parameters for M random effects
        M_re_model: (n_stocks x n_regions) 1: no RE, 2: iid/ar1 RE for age (constant by year), 3: iid/ar1 RE year (constant by age), 4: 2D iid/AR1 (age, year) (not yet 5: 3DAR1 (age,year,region))
           M_model: 1: M= f(age), 2: M = f(WAA)
              M_re: random effects for 
            n_M_re: (n_stocks x n_regions) number of random effects (e.g., number of age classes or number of M=f(WAA) pars) each year), max(n_M_re)<=n_ages)
         years_use: is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = M_re.dim(0);
  int n_regions = M_re.dim(1);
  //int n_ages = M_re.dim(3);
  int n_y = years_use.size();
  //int n_y = M_re.dim(2);
  array<Type> sim_M_re = M_re;
  //array<Type> sim_M_re(n_stocks,n_regions,n_y,n_ages);// = M_re;
  //sim_M_re.setZero();
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r< n_regions; r++){
    if(M_re_model(s,r) > 1){ // random effects on M, M_re = 2D AR1 deviations on M(year,age), dim = n_years x n_M_re(s,r)
      Type sigma_M = exp(M_repars(s,r,0));
      Type rho_M_a = geninvlogit(M_repars(s,r,1),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy
      Type rho_M_y = geninvlogit(M_repars(s,r,2),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy
      Type Sigma_M;
      // likelihood of M deviations, M_re
      array<Type> M_re_r_s(n_y,n_M_re(s,r));
      M_re_r_s.setZero();
      if(M_re_model(s,r) == 4){ //2D AR1: age, year (could be iid)
        Sigma_M = pow(pow(sigma_M,2) / ((1-pow(rho_M_y,2))*(1-pow(rho_M_a,2))),0.5);
        SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)).simulate(M_re_r_s); // must be array, not matrix!
        for(int y = 0; y < n_y; y++)for(int a = 0; a < n_M_re(s,r); a++) sim_M_re(s,r,years_use(y),a) = M_re_r_s(y,a)*Sigma_M;
      } else {
        if(M_re_model(s,r) == 2){ // 1D ar1_a
          vector<Type> Mre0 = M_re_r_s.matrix().row(0);
          Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_a,2)),0.5);
          AR1(rho_M_a).simulate(Mre0);
          //all years mapped to the same age-specific RE
          for(int y = 0; y < n_y; y++) for(int a = 0; a < n_M_re(s,r); a++) sim_M_re(s,r,years_use(y),a) = Mre0(a)*Sigma_M;
        } else { // M_re_model(s,r) = 3, 1D ar1_y
          vector<Type> Mre0 = M_re_r_s.matrix().col(0);
          Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_y,2)),0.5);
          AR1(rho_M_a).simulate(Mre0);
          //all ages mapped to the same annual RE
          for(int y = 0; y < n_y; y++) for(int a = 0; a < n_M_re(s,r); a++) sim_M_re(s,r,years_use(y),a) = Mre0(y)*Sigma_M;
        }
      }
    }
  }
  return(sim_M_re);
}
//done

//nll for prior on M(WAA) (log) exponent
template<class Type>
matrix<Type> get_nll_log_b(int log_b_model, matrix<Type>log_b, int bias_correct){
  /* 
    get any nll components for prior on log_b of M_model = 2 where M = a * W ^b
       log_b_model: model for posterior log_b. 1: stock and region constant, 2: region-constant, 3: stock constant, 4: differ by stock and regio
             log_b: random effect parameters for posterior log_b for each stock,region
      bias_correct: 0/1 whether to bias-correct log(b) prior mean 
  */
  int n_stocks = log_b.rows();
  int n_regions = log_b.cols();
  matrix<Type> nll(n_stocks,n_regions);
  nll.setZero();
  Type mu = log(0.305);
  if(bias_correct) mu -= 0.5*exp(2*log(0.08));
  if(log_b_model == 1){ //constant across stocks, regions
    nll(0,0) = -dnorm(log_b(0,0), mu, Type(0.08),1);
    //for(int r = 0; r < n_regions; r++) for(int s = 0; s < n_stocks; s++) sim_log_b(r,s) = sim;
  }
  if(log_b_model == 2){ //constant across regions, differ by stock
    for(int r = 0; r < n_regions; r++) {
      nll(0,r) = -dnorm(log_b(0,r), mu, Type(0.08), 1);
      //for(int s = 0; s < n_stocks; s++) sim_log_b(r,s) = sim;
    }
  }
  if(log_b_model == 3){ //constant across stocks, differ by region
    for(int s = 0; s < n_stocks; s++) {
      nll(s,0) = -dnorm(log_b(s,0), mu, Type(0.08), 1);
      //for(int r = 0; r < n_regions; r++) sim_log_b(r,s) = sim;
    }
  }
  if(log_b_model == 4){ //differ by stocks and region
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
      nll(s,r) = -dnorm(log_b(s,r), mu, Type(0.08), 1);
      //sim_log_b(r,s) = rnorm(mu, Type(0.08));
    }
  }
  return(nll);
}
//done

//simulate log exponent for M(WAA)
template<class Type>
matrix<Type> simulate_log_b(int log_b_model, matrix<Type>log_b, int bias_correct){
  /* 
    simulate any RE for log_b of M_model = 2 where M = a * W ^b
       log_b_model: model for posterior log_b. 1: stock and region constant, 2: region-constant, 3: stock constant, 4: differ by stock and regio
             log_b: random effect parameters for posterior log_b for each stock,region
      bias_correct: 0/1 whether to bias-correct log(b) prior mean 
  */
  int n_stocks = log_b.rows();
  int n_regions = log_b.cols();
  matrix<Type> sim_log_b = log_b;
  Type mu = log(0.305);
  if(bias_correct == 1) mu -= 0.5*exp(2*log(0.08));
  if(log_b_model == 1){ //constant across stocks, regions
    Type sim = rnorm(mu, Type(0.08));
    for(int r = 0; r < n_regions; r++) for(int s = 0; s < n_stocks; s++) sim_log_b(r,s) = sim;
  }
  if(log_b_model == 2){ //constant across regions, differ by stock
    for(int r = 0; r < n_regions; r++) {
      Type sim = rnorm(mu, Type(0.08));
      for(int s = 0; s < n_stocks; s++) sim_log_b(r,s) = sim;
    }
  }
  if(log_b_model == 3){ //constant across stocks, differ by region
    for(int s = 0; s < n_stocks; s++) {
      Type sim = rnorm(mu, Type(0.08));
      for(int r = 0; r < n_regions; r++) sim_log_b(r,s) = sim;
    }
  }
  if(log_b_model == 4){ //differ by stocks and region
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
      sim_log_b(r,s) = rnorm(mu, Type(0.08));
    }
  }
  return(sim_log_b);
}
//done

template<class Type>
array<Type> get_log_M(array<Type>M_re, array<int> M_re_index, int M_model, int n_years_model, array<Type> Mpars, matrix<Type> log_b, array<Type> waa, 
  vector<int> waa_pointer, array<Type> Ecov_lm, array<int> Ecov_how, int n_years_proj, int proj_M_opt, vector<int> avg_years_ind){
  /* 
    provides log_M (base) components that are density-independent.
    (log) mortality-at-age (MAA)
               M_re: n_stocks x n_regions x n_years x n_ages
         M_re_index: (n_stocks x n_regions x n_ages): indicators (consistent with n_M_re) indicating which ages to apply M_re to.
            M_model: 1: M= f(age), 2: M = f(WAA)
      n_years_model: number of years in the population model (needed to tell when to start projection options)
                Mpars: (n_stocks x n_regions x n_ages) mean M-at-age, fixed effects, if M_model = 1 (age-constant, age-specific) all ages are used. If M_model = 2 (M = a * W ^b) just first one used 
              log_b: random effect parameters for posterior log_b (M_model = 2 M = a * W ^b) for each stock,region
                waa: (at least(n_fleets + n_indices + n_stocks + 1(totcatch)) x n_years x n_ages_model) weight at age
        waa_pointer: n_stocks, which index in first dimension to use if M_model = 2
            Ecov_lm: (n_stocks, n_regions, n_years_pop, n_Ecov) linear predictor for any Ecov effects on log_M
           Ecov_how: n_Ecov x n_stocks x n_ages x n_regions: 0/1 values indicating to use effects on natural mortality at age.
       n_years_proj: number of projection years
         proj_M_opt: 1 = continue M_re (check for time-varying M_re on R side), 2 = average M (over avg_years_ind)
      avg_years_ind: model year indices (TMB, starts @ 0) to use for averaging MAA, waa, maturity, and F (if use.avgF = TRUE)
  */
  int n_stocks = M_re.dim(0);
  int n_regions = M_re.dim(1);
  int n_ages = M_re.dim(3);
  int n_years = M_re.dim(2);
  array<Type> log_M(n_stocks,n_regions,n_years,n_ages);
  log_M.setZero();
  for(int r = 0; r< n_regions; r++) for(int s = 0; s< n_stocks; s++){    
    // Mpars should be mapped properly
    for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) {
      log_M(s,r,y,a) = Mpars(s,r,a);
      if(M_model == 1){ // age-specific or constant M (maybe by stock)
        if(M_re_index(s,r,a)>0) log_M(s,r,y,a) += M_re(s,r,y,M_re_index(s,r,a)-1);
      } else { // M_model = 2, M is allometric function of weight
        log_M(s,r,y,a) += M_re(s,r,y,0) - exp(log_b(s,r) + M_re(s,r,y,1)) * log(waa(waa_pointer(s)-1,y,a));
      }
    } 
    // add ecov effect on M
    for(int i=0; i < Ecov_how.dim(0); i++) for(int a = 0; a < n_ages; a++) {
      if(Ecov_how(i,s,a,r)==1) for(int y = 0; y < n_years_model; y++) log_M(s,r,y,a) += Ecov_lm(s,r,a,y,i);
    }
    
    // add to MAA in projection years
    if(n_years_proj>0){ 
      int n_toavg = avg_years_ind.size();
      if(proj_M_opt == 2){ // use average MAA over avg.yrs 
        matrix<Type> MAA_toavg(n_toavg,n_ages);
        for(int a = 0; a < n_ages; a++){
          for(int i = 0; i < n_toavg; i++){
            MAA_toavg(i,a) = exp(log_M(s,r,avg_years_ind(i),a));
          }
        }
        vector<Type> MAA_proj = MAA_toavg.colwise().mean();
        for(int y = n_years_model; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
          log_M(s,r,y,a) = log(MAA_proj(a));
          //MAA.row(y) = MAA_proj;
        }
      } else { // proj_M_opt == 1, use M_re and/or ecov_lm in projection years
        for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < n_years; y++) {
          log_M(s,r,y,a) = Mpars(s,r,a);
          if((M_model == 1)){ // age-specific or constant M (maybe by stock)
            if(M_re_index(s,r,a)>0) {
              log_M(s,r,y,a) += M_re(s,r,y,M_re_index(s,r,a)-1);
            }
          } else { // M_model = 2, M is allometric function of weight
            for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < n_years; y++) {
              log_M(s,r,y,a) += M_re(s,r,y,0) - exp(log_b(s,r) + M_re(s,r,y,1)) * log(waa(waa_pointer(s)-1,y,a));
            }
          }
        }
      }
      // add ecov effect on M
      for(int i=0; i < Ecov_how.dim(0); i++) for(int a = 0; a < n_ages; a++) {
        if(Ecov_how(i,s,a,r)==1) for(int y = n_years_model; y < n_years; y++) log_M(s,r,y,a) += Ecov_lm(s,r,a,y,i);
      }
    }
  }
  return log_M; 
}

template<class Type>
array<Type> get_log_avg_M(array<Type> log_M, vector<int> years_M) {
  array<Type> M(log_M.dim(0),log_M.dim(1),log_M.dim(3));
  array<Type> log_avg_M(log_M.dim(0),log_M.dim(1),log_M.dim(3));
  for(int s = 0; s < log_M.dim(0); s++) for(int r = 0; r < log_M.dim(1); r++) for(int a = 0; a < log_M.dim(3); a++) {
    for(int y = 0; y < years_M.size(); y++) M(s,r,a) += exp(log_M(s,r,years_M(y),a))/Type(years_M.size());
    log_avg_M(s,r,a) = log(M(s,r,a)); 
  }
  return log_avg_M;
}

//extract array of log_M parameters for a given year
template <class Type>
array<Type> get_log_M_y(int y, array<Type> log_M, int do_log = 1){
  int n_stocks = log_M.dim(0);
  int n_ages = log_M.dim(3);
  int n_regions = log_M.dim(1);
  array<Type> log_M_y(n_stocks,n_regions,n_ages);
  log_M_y.setZero();

  for(int s = 0; s< n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++){
    if(do_log) {
      log_M_y(s,r,a) = log_M(s,r,y,a);
    } else {
      log_M_y(s,r,a) = exp(log_M(s,r,y,a));
    }
  }
  return log_M_y;
}

//exponentiate log_M
template <class Type>
array<Type> get_MAA(array<Type> log_M){
  array<Type> MAA(log_M.dim(0),log_M.dim(1),log_M.dim(2),log_M.dim(3));
  MAA.setZero();

  for(int s = 0; s< log_M.dim(0); s++) for(int r = 0; r < log_M.dim(1); r++) {
    for(int y = 0; y < log_M.dim(2); y++) for(int a = 0; a < log_M.dim(3); a++){
      MAA(s,r,y,a) = exp(log_M(s,r,y,a));
    }
  }
  return MAA;
}
