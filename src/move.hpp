// mu_model: 
// 1 = constant across stocks, ages, time (n_seasons x n_regions x (n_regions -1) pars). 
// 2 = differ by age (n_seasons x n_regions x (n_regions -1) fixed effects, n_ages random effects for each). 
// 3 = differ by year (n_seasons x n_regions x (n_regions -1) fixed effects, n_years random effects for each)
// 4 = differ by age,year (n_seasons x n_regions x (n_regions -1) fixed effects, n_years x n_ages random effects for each)
// 5 = differ by stock (n_seasons x n_stocks x n_regions x (n_regions -1) pars). 
// 6 = differ by stock, age (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_ages random effects for each). 
// 7 = differ by stock, year (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_years random effects for each)
// 8 = differ by stock, age,year (n_seasons x n_stocks x n_regions x (n_regions -1) fixed effects, n_years x n_ages random effects for each)

template <class Type>
array<Type> get_nll_mu_prior(array<Type> mu_prior_re, array<Type> trans_mu, array<Type> trans_mu_prior_sigma, array<int> use_mu_pior, int mu_model){
  int n_stocks = trans_mu.dim(0);
  int n_ages = trans_mu.dim(1);
  int n_seasons = trans_mu.dim(2);
  int n_regions = trans_mu.dim(3);
  array<Type> nll(n_stocks,n_ages, n_seasons, n_regions, n_regions-1);
  nll.setZero();
  if((mu_model == 1) || (mu_model == 2) || (mu_model == 3) || (mu_model == 4) ){
    for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(0,0,t,r,rr)==1) nll(0,0,t,r,rr) -= dnorm(mu_prior_re(0,0,t,r,rr), trans_mu(0,0,t,r,rr), trans_mu_prior_sigma(0,0,t,r,rr),1);
    }
  }
  if((mu_model == 5) || (mu_model == 6) || (mu_model == 7) || (mu_model == 8)){
    for(int s = 0; s < n_stocks; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(s,0,t,r,rr)==1) nll(s,0,t,r,rr) -= dnorm(mu_prior_re(s,0,t,r,rr), trans_mu(s,0,t,r,rr), trans_mu_prior_sigma(s,0,t,r,rr),1);
    }
  }
  return(nll);
}

template <class Type>
array<Type> simulate_mu_prior_re(array<Type> mu_prior_re, array<Type> trans_mu, array<Type> trans_mu_prior_sigma, array<int> use_mu_prior, int mu_model){
  int n_stocks = trans_mu.dim(0);
  int n_ages = trans_mu.dim(1);
  int n_seasons = trans_mu.dim(2);
  int n_regions = trans_mu.dim(3);
  array<Type> sim_mu_prior_re(n_stocks,n_ages, n_seasons, n_regions, n_regions-1);
  sim_mu_prior_re.setZero();
  if((mu_model == 1) || (mu_model == 2) || (mu_model == 3) || (mu_model == 4) ){
    for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(0,0,t,r,rr)==1) {
        sim_mu_prior_re(0,0,t,r,rr) = rnorm(trans_mu(0,0,t,r,rr), trans_mu_prior_sigma(0,0,t,r,rr));
      }
    }
  }
  if((mu_model == 5) || (mu_model == 6) || (mu_model == 7) || (mu_model == 8)){
    for(int s = 0; s < n_stocks; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(s,0,t,r,rr)==1){
        sim_mu_prior_re(s,0,t,r,rr) = rnorm(trans_mu(s,0,t,r,rr), trans_mu_prior_sigma(s,0,t,r,rr));
      }
    }
  }
  return(sim_mu_prior_re);
}

template <class Type>
array<Type> get_nll_mu_re(array<Type> mu_repars, array<Type> mu_re, int mu_model){
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_y = mu_re.dim(3);
  int n_regions = mu_re.dim(4);
  array<Type> nll(n_stocks,n_age,n_seasons,n_regions,n_regions-1);
  nll.setZero();
  if((mu_model != 1) && (mu_model != 5)){ //some type of random effects
    int ns = 1;
    if(mu_model > 5) ns = n_stocks; 
    if((mu_model == 2) || (mu_model == 6)){ // age random effects
      for(int s = 0; s < ns; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++){
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = rho_trans(mu_repars(s,0,t,r,1));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_a,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_a(n_ages);
          for(int a = 0; a < n_ages; a++) mu_re_a(a) = mu_re(s,a,t,0,r,rr);
          nll(s,0,t,r,rr) += SCALE(AR1(rho_mu_a), Sigma_MU)(mu_re_a);
        }
      }
    }
    if((mu_model == 3) || (mu_model == 7)){ // year random effects
      for(int s = 0; s < ns; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_y = rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_y,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_y(n_y);
          for(int y = 0; y< n_y; y++) mu_re_y(y) = mu_re(s,0,t,y,r,rr);
          nll(s,0,t,r,rr) += SCALE(AR1(rho_mu_y), Sigma_MU)(mu_re_y);
        }
      }
    }
    if((mu_model == 4) || (mu_model == 8)){ // age and year random effects
      for(int s = 0; s < ns; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = rho_trans(mu_repars(s,0,t,r,1));
        Type rho_mu_y = rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_y,2)),0.5); //marginal sd
        Sigma_MU = pow(pow(sigma_mu,2) / ((1-pow(rho_mu_y,2))*(1-pow(rho_mu_a,2))),0.5);//marginal variance
        for(int rr = 0; rr < n_regions-1; rr++){
          array<Type> mu_re_ya(n_y,n_ages);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) mu_re_ya(y,a) = mu_re(s,a,t,y,r,rr);
          nll(s,0,t,r,rr) += SCALE(SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)), Sigma_M)(mu_re_ya); // must be array, not matrix!
        }
      }
    }
  }
  return(nll);
}

template <class Type>
array<Type> simulate_mu_re(array<Type> mu_repars, array<Type> mu_re, int mu_model){
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_y = mu_re.dim(3);
  int n_regions = mu_re.dim(4);
  array<Type> sim_mu_re(n_stocks,n_age,n_seasons,n_regions,n_regions-1);
  sim_mu_re.setZero();
  if((mu_model != 1) && (mu_model != 5)){ //some type of random effects
    int ns = 1;
    if(mu_model > 5) ns = n_stocks; 
    if((mu_model == 2) || (mu_model == 6)){ // age random effects
      for(int s = 0; s < ns; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++){
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = rho_trans(mu_repars(s,0,t,r,1));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_a,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_a(n_ages);
          AR1(rho_mu_a).simulate(mu_re_a);
          for(int a = 0; a < n_ages; a++) sim_mu_re(s,a,t,0,r,rr) = mu_re_a(a) * Sigma_MU;
        }
      }
    }
    if((mu_model == 3) || (mu_model == 7)){ // year random effects
      for(int s = 0; s < ns; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_y = rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_y,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_y(n_y);
          AR1(rho_mu_y).simulate(mu_re_y);
          for(int y = 0; y< n_y; y++) sim_mu_re(s,0,t,y,r,rr) = mu_re_y(y) * Sigma_MU;
        }
      }
    }
    if((mu_model == 4) || (mu_model == 8)){ // age and year random effects
      for(int s = 0; s < ns; s++) for(int t = 0, t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = rho_trans(mu_repars(s,0,t,r,1));
        Type rho_mu_y = rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_y,2)),0.5); //marginal sd
        Sigma_MU = pow(pow(sigma_mu,2) / ((1-pow(rho_mu_y,2))*(1-pow(rho_mu_a,2))),0.5);//marginal variance
        for(int rr = 0; rr < n_regions-1; rr++){
          array<Type> mu_re_ya(n_y,n_ages);
          SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)).simulate(mu_re_ya);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++)  sim_mu_re(s,a,t,y,r,rr) = mu_re_ya(y,a) * Sigma_MU;
        }
      }
    }
  }
  return(sim_mu_re);
}


template <class Type>
additive_ln_transform(vector<Type> x, int region, vector<int> do_move, int must_move){
  //use additive transformation (e.g., logistic-normal model)
  //ensures that probabilities of moving and staying add to 1
  int D = x.size()+1;
  vector<Type> y(D);
  y.setZero();
  int j = 0;
  for(int i = 0; i < D; i++) {
    if(i != region) {
      if(do_move(i)==1) y(i) = exp(x(j)); //else prob of moving will be 0.
      j++;
    } 
    if(i ==region) { //prob of staying will be 1- prob of moving
      if(must_move==0) {
        y(i) = 1.0;
      } //else y(i) = 0.0 already. However, must make sure that trans_mu is fixed at 0 for one of the do_move regions to make returned y sum to 1.
    }
  }
  y /= sum(y);
  return(y);
}

//provides transformed mu (good for sdreporting)
template<class Type>
array<Type> get_trans_mu_base(array<Type> trans_mu, array<Type>mu_re, array<Type> mu_prior_re, array<int> use_mu_prior, int mu_model,
  array<Type> Ecov_lm, array<int> use_Ecov){
  // Construct base mu-at-age 
  // currently continues random processes in any projection years!
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_regions = mu_re.dim(4);
  int ny = mu_re.dim(3);
  //array<Type> Ecov_lm_mu(n_stocks, n_regions-1, n_ages, n_seasons, n_years_model + n_years_proj, n_Ecov);
  array<Type> trans_mu_base(n_stocks,n_ages,n_seasons,ny,n_regions, n_regions-1);
  trans_mu_base.setZero();
  for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0, t < n_seasons; t++) for(int y = 0; y < n_y; y++){
    for(int r = 0; r< n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++) {
      if(use_mu_prior(s,a,t,r,rr)==1) {
        trans_mu_base(s,a,t,y,r,rr) += mu_prior_re(s,a,t,r,rr);
      } else {
        trans_mu_base(s,a,t,y,r,rr) += trans_mu(s,a,t,r,rr);
      }
      if((mu_model != 1) && (mu_model != 5)){ //some type of random effects
        if(mu_model == 2) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,a,t,0,r,rr); // age random effects
        if(mu_model == 6) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,a,t,0,r,rr);// age random effects by stock
        if(mu_model == 3) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,0,t,y,r,rr); // year random effects
        if(mu_model == 7) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,0,t,y,r,rr);// year random effects by stock
        if(mu_model == 4) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,a,t,y,r,rr); // age,year random effects
        if(mu_model == 8) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,a,t,y,r,rr);// age,year random effects by stock
      }
      for(int i=0; i < use_Ecov.dim(0); i++) if(use_Ecov(i,s,a,t,r,rr) == 1) trans_mu_base(s,a,t,y,r,rr) += Ecov_lm(s,a,t,y,r,rr,i); //will be 0 if not used
    }
  }
  //no projections options for mu. Just forecast any random or Ecov effects. Otherwise constant mu is the same as during model period.
  return(trans_mu_base); 
}

//n_regions x n_regions movement matrix
template <class Type>
matrix<Type> get_mu_matrix(int stock, int age, int season, int year, vector<int> mig_type, array<int> can_move, array<int> must_move, array<Type> trans_mu_base){
  int n_regions = trans_mu_base.dim(4);
  matrix<Type> mu(n_regions,n_regions);
  mu.setZero();
  if(mig_type(stock) == 0) //migration is instantaneous after survival and mortality, so P is easy.
  {
    for(int r = 0; r < n_regions; r++) {
      vector<Type> trans_par(n_regions-1);
      for(int rr = 0; rr < n_regions-1; rr++) trans_par(rr) = trans_mu_base(stock,age,season,year,r,rr);
      for(int j = 0; j < n_regions; j++) do_move(j) = can_move(s,a,t,r,j);
      vector<Type> pmove = additive_transform(trans_par, r, do_move, must_move(s,a,t,r));
      for(int j = 0; j < n_regions; j++) mu(r,j) = pmove(j);
    }
  }
  if(mig_type(stock) == 1) //migration occurs continuously during interval, return infinitesimal generator.
  {
    for(int r = 0; r < n_regions; r++) {
      int k = 0;
      for(int j = 0; j < n_regions; j++){ 
        if(j!=i) {
          k++; //max k = n_regions -1 (-1)
          if(can_move(i,j)==1) mu(i,j) = exp(trans_mu_base(stock,age,season,year,r,k)); //log of transition intensities
        }
      }
    }
    for(int r = 0; r< n_regions; r++) mu(r,r) = -(mu.row(r)).sum(); //hazard
  }
  return(mu);
}

//all movement matrices
template <class Type>
array<Type> get_mu(array<Type> trans_mu_base, array<int> can_move, 
  array<int> must_move, vector<int> mig_type){
  //additive logistic-normal transform for probability of movement out of regions only; prob staying is 1 - sum(prob move)
  int n_stocks = trans_mu_base.dim(0);
  int n_ages = trans_mu_base.dim(1);
  int n_seasons = trans_mu_base.dim(2);
  int n_y = trans_mu_base.dim(3);
  int n_regions = trans_mu_base.dim(4);
  array<Type> mu(m_stocks,n_ages, n_seasons, n_years, n_regions,n_regions);
  mu.setZero();
  for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0, t < n_seasons; t++) for(int y = 0; y < n_y; y++){
    matrix<Type> mu_y = get_mu_matrix(s,a,t,y,mig_type,can_move,must_move,trans_mu_base);
    for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) mu(s,a,t,y,r,rr) = mu_y(r,rr);
  }
  return(mu);
}
  
