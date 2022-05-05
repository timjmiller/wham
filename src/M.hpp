template<class Type>
matrix<Type> get_nll_M(array<Type> M_repars, vector<int> M_re_model, vector<int> M_model, array<Type>M_re, matrix<int> n_M_re){
  int n_stocks = M_re.dim(0);
  int n_regions = M_re.dim(1);
  int n_ages = M_re.dim(3);
  int n_y = M_re.dim(2);
  matrix<Type> nll_M(n_stocks,n_regions);
  nll_M.setZero();
  for(int r = 0; r< n_regions; r++){
    if(M_re_model(r) > 1) // random effects on M, M_re = 2D AR1 deviations on M(year,age), dim = n_years x n_M_re(s,r)
    {
      int n_stock_re = 1;
      if(M_model(r)>3) n_stock_re = n_stocks; //different M models for each stock in each region.
      for(int s = 0; s < n_stock_re; s++){
        Type sigma_M = exp(M_repars(s,r,0));
        Type rho_M_a = rho_trans(M_repars(s,r,1));
        Type rho_M_y = rho_trans(M_repars(s,r,2));
        Type Sigma_M;
        // likelihood of M deviations, M_re
        array<Type> M_re_r_s(n_y,n_M_re(s,r));
        M_re_r_s.setZero();
        for(int y = 0; y < n_y; y++)for(int a = 0; a < n_M_re(s,r); a++) M_re_r_s(y,a) = M_re(s,r,y,a); //first n_M_re(s,r) columns
        if((M_re_model(r) == 2) | (M_re_model(r) == 5)){ //2D AR1: age, year
          Sigma_M = pow(pow(sigma_M,2) / ((1-pow(rho_M_y,2))*(1-pow(rho_M_a,2))),0.5);
          nll_M += SCALE(SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)), Sigma_M)(M_re_r_s); // must be array, not matrix!
        } else {
          if(M_re_model == 3){ // 1D ar1_a
            vector<Type> Mre0 = M_re_r_s.row(0);
            Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_a,2)),0.5);
            nll_M += SCALE(AR1(rho_M_a), Sigma_M)(Mre0);
          } else { // M_re_model = 4, 1D ar1_y
            vector<Type> Mre0 = M_re_r_s.col(0); //just first column
            Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_y,2)),0.5);
            nll_M += SCALE(AR1(rho_M_y), Sigma_M)(Mre0);
          }
        }
      }
    }
  }
  return(nll_M);
}

template<class Type>
array<Type> simulate_M_re(array<Type> M_repars, vector<int> M_re_model, vector<int> M_model, array<Type>M_re, matrix<int> n_M_re){
  int n_stocks = M_re.dim(0);
  int n_regions = M_re.dim(1);
  /int n_ages = M_re.dim(3);
  int ny = M_re.dim(2);
  array<Type> sim_M_re = M_re;
  for(int r = 0; r< n_regions; r++){
    if(M_re_model(r) > 1) // random effects on M, M_re = 2D AR1 deviations on M(year,age), dim = n_years x n_M_re(s,r)
    {
      int n_stock_re = 1;
      if(M_model(r)>3) n_stock_re = n_stocks; //different M models for each stock in each region.
      for(int s = 0; s < n_stock_re; s++){
        Type sigma_M = exp(M_repars(s,r,0));
        Type rho_M_a = rho_trans(M_repars(s,r,1));
        Type rho_M_y = rho_trans(M_repars(s,r,2));
        Type Sigma_M;
        // likelihood of M deviations, M_re
        array<Type> M_re_r_s(n_y,n_M_re(s,r));
        M_re_r_s.setZero();
        if((M_re_model(r) == 2) | (M_re_model(r) == 5)){ //2D AR1: age, year
          Sigma_M = pow(pow(sigma_M,2) / ((1-pow(rho_M_y,2))*(1-pow(rho_M_a,2))),0.5);
          SEPARABLE(AR1(rho_M_a),AR1(rho_M_y)).simulate(M_re_r_s); // must be array, not matrix!
          for(int y = 0; y < n_y; y++)for(int a = 0; a < n_M_re(s,r); a++) sim_M_re(s,r,y,a) = M_re_r_s(y,a)*Sigma_M;
        } else {
          if(M_re_model == 3){ // 1D ar1_a
            vector<Type> Mre0 = M_re_r_s.row(0);
            Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_a,2)),0.5);
            AR1(rho_M_a).simulate(Mre0);
            for(int a = 0; a < n_M_re(s,r); a++) sim_M_re(s,r,0,a) = Mre0(a)*Sigma_M;
          } else { // M_re_model = 4, 1D ar1_y
            vector<Type> Mre0 = M_re_r_s.col(0);
            Sigma_M = pow(pow(sigma_M,2) / (1-pow(rho_M_y,2)),0.5);
            AR1(rho_M_a).simulate(Mre0);
            for(int y = 0; y < n_y; y++) sim_M_re(s,r,y,0) = Mre0(y)*Sigma_M;
          }
        }
      }
    }
  }
  return(sim_M_re);
}

//provides log_M components that are density-independent.
template<class Type>
array<Type> get_log_M_base(array<Type>M_re, int M_model, int n_years_model, array<Type> M_a, Type log_b, array<Type> waa, 
  int waa_pointer, array<Type> Ecov_lm, array<int> use_Ecov, int do_proj, int proj_M_opt, vector<int> avg_years_ind){
  // Construct base (log) mortality-at-age (MAA)
  int n_stocks = M_re.dim(0);
  int n_regions = M_re.dim(1);
  int n_ages = M_re.dim(3);
  int ny = M_re.dim(2);
  array<Type> log_M_base(n_stocks,n_regions,ny,n_ages);
  log_M_base.setZero();
  for(int r = 0; r< n_regions; r++) for(int s = 0; s< n_stocks; s++){
    if(M_model == 2){ // age-specific M
      for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) log_M_base(s,r,y,a) = M_a(s,r,a) + M_re(s,r,y,a);   
    } else {
      if(M_model == 1){ // constant M
        for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) log_M_base(s,r,y,a) = M_a(s,r,0) + M_re(s,r,y,a);
      } else { // M_model = 3, M is allometric function of weight
        for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_years_model; y++) log_M_base(s,r,y,a) = M_a(s,r,0) + M_re(s,r,y,a) - exp(log_b) * log(waa(waa_pointer-1,y,a));
      }
    }
    // add ecov effect on M (by year, shared across ages)
    for(int i=0; i < use_Ecov.dim(0); i++) for(int a = 0; a < n_ages; a++) {
      if(use_Ecov(i,s,a,r) == 1) for(int y = 0; y < n_years_model; y++) log_M_base(s,r,y,a) += Ecov_lm(s,r,a,y,i);
    }
    
    // add to MAA in projection years
    if(do_proj == 1){ 
      int n_toavg = avg_years_ind.size();
      if(proj_M_opt == 2){ // use average MAA over avg.yrs 
        matrix<Type> MAA_toavg(n_toavg,n_ages);
        for(int a = 0; a < n_ages; a++){
          for(int i = 0; i < n_toavg; i++){
            MAA_toavg(i,a) = exp(log_M_base(s,r,avg_years_ind(i),a));
          }
        }
        vector<Type> MAA_proj = MAA_toavg.colwise().mean();
        for(int y = n_years_model; y < ny; y++) for(int a = 0; a < n_ages; a++) {
          log_M_base(s,r,y,a) = log(MAA_proj(a));
          //MAA.row(y) = MAA_proj;
        }
      } else { // proj_M_opt == 1, use M_re and/or ecov_lm in projection years
        if(M_model == 2){ // age-specific M
          for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < ny; y++) log_M_base(s,r,y,a) = M_a(s,r,a) + M_re(s,r,y,a);   
        } else {
          if(M_model == 1){ // constant M
            for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < ny; y++) log_M_base(s,r,y,a) = M_a(s,r,0) + M_re(s,r,y,a);
          } else { // M_model = 3, M is allometric function of weight
            for(int a = 0; a < n_ages; a++) for(int y = n_years_model; y < ny; y++) log_M_base(s,r,y,a) = M_a(s,r,0) + M_re(s,r,y,a) - exp(log_b) * log(waa(waa_pointer-1,y,a));
          }
        }
      }
    }
  }
  return(log_M_base); 
}

