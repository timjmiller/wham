// mu_model(r,rr): n_regions x (n_regions -1) 
// 1 = constant across stocks, ages, time (1 fixed effect for r,rr). 
// 2 = differ by age (1 fixed effect, n_ages random effects for r,rr). 
// 3 = differ by year, (1 fixed effect, n_years random effects for r,rr)
// 4 = differ by age,year (1 fixed effect, n_years x n_ages random effects for r,rr)
// 5 = differ by stock (n_stocks fixed effects for r,rr). 
// 6 = differ by stock, age (n_stocks fixed effects, n_ages random effects for r,rr). 
// 7 = differ by stock, year (n_stocks fixed effects, n_years random effects for r,rr)
// 8 = differ by stock, age,year (n_stocks fixed effects, n_years x n_ages random effects for r,rr)
// 9 = differ by season (n_seasons fixed effects for r,rr). 
// 10 = differ by season,age (n_seasons fixed effects, n_ages random effects for r,rr). 
// 11 = differ by season,year (n_seasons fixed effects, n_years random effects for r,rr)
// 12 = differ by season,age,year (n_seasons fixed effects, n_years x n_ages random effects for r,rr)
// 13 = differ by stock, season (n_stocks x n_seasons fixed effects for r,rr). 
// 14 = differ by stock, season, age (n_stocks x n_seasons fixed effects, n_ages random effects for r,rr). 
// 15 = differ by stock, season, year (n_stocks x n_seasons fixed effects, n_years random effects for r,rr)
// 16 = differ by stock, season, age,year (n_stocks x n_seasons fixed effects, n_years x n_ages random effects for r,rr)
template <class Type>
array<Type> get_nll_mu_prior(array<Type> mu_prior_re, array<Type> trans_mu, array<Type> trans_mu_prior_sigma, array<int> use_mu_prior, matrix<int> mu_model){
  /* 
    get any nll components for priors/posteriors on mu parameters.
               mu_prior_re: n_stocks x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
                  trans_mu: n_stocks x n_seasons x n_regions x n_regions-1 (mean) movement parameters
      trans_mu_prior_sigma: n_stocks x n_seasons x n_regions x n_regions-1 sd of prior on transformed movement parameters
              use_mu_prior: n_stocks x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
                  mu_model: n_regions x n_regions-1. see definitions at top of move.hpp.
  */
  int n_stocks = trans_mu.dim(0);
  //int n_ages = trans_mu.dim(1);
  int n_seasons = trans_mu.dim(1);
  int n_regions = trans_mu.dim(2);
  array<Type> nll(n_stocks,n_seasons, n_regions, n_regions-1);
  nll.setZero();
  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
    if((mu_model(r,rr) > 0) & (mu_model(r,rr) <= 4)){ //constant
      if(use_mu_prior(0,0,r,rr)) nll(0,0,r,rr) -= dnorm(mu_prior_re(0,0,r,rr), trans_mu(0,0,r,rr), trans_mu_prior_sigma(0,0,r,rr),1);
    }
    if((mu_model(r,rr) > 4) & (mu_model(r,rr) <= 8)){ //stock
      for(int s = 0; s < n_stocks; s++){
        if(use_mu_prior(s,0,r,rr)) nll(s,0,r,rr) -= dnorm(mu_prior_re(s,0,r,rr), trans_mu(s,0,r,rr), trans_mu_prior_sigma(s,0,r,rr),1);
      }
    }
    if((mu_model(r,rr) > 8) & (mu_model(r,rr) <= 12)){ //season
      for(int t = 0; t < n_seasons; t++){
        if(use_mu_prior(0,t,r,rr)) nll(0,t,r,rr) -= dnorm(mu_prior_re(0,t,r,rr), trans_mu(0,t,r,rr), trans_mu_prior_sigma(0,t,r,rr),1);
      }
    }
    if((mu_model(r,rr) > 12) & (mu_model(r,rr) <= 16)){ //stock,season
      for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++){
        if(use_mu_prior(s,t,r,rr)) nll(s,t,r,rr) -= dnorm(mu_prior_re(s,t,r,rr), trans_mu(s,t,r,rr), trans_mu_prior_sigma(s,t,r,rr),1);
      }
    }
  }
  return(nll);
}
//done

template <class Type>
array<Type> simulate_mu_prior_re(array<Type> mu_prior_re, array<Type> trans_mu, array<Type> trans_mu_prior_sigma, array<int> use_mu_prior, matrix<int> mu_model){
  /* 
    simulate and RE for priors on mu parameters.
               mu_prior_re: n_stocks x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
                  trans_mu: n_stocks x n_seasons x n_regions x n_regions-1 (mean) movement parameters
      trans_mu_prior_sigma: n_stocks x n_seasons x n_regions x n_regions-1 (mean) movement parameters
              use_mu_prior: n_stocks x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
                  mu_model: n_regions x n_regions-1. see definitions at top of move.hpp.
  */
  int n_stocks = trans_mu.dim(0);
  //int n_ages = trans_mu.dim(1);
  int n_seasons = trans_mu.dim(1);
  int n_regions = trans_mu.dim(2);
  array<Type> sim_mu_prior_re(n_stocks, n_seasons, n_regions, n_regions-1);
  sim_mu_prior_re.setZero();
    array<Type> sims(n_regions,n_regions-1);
    sims.setZero();
  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
    if((mu_model(r,rr) > 0) & (mu_model(r,rr) <= 4)) {//constant
      if(use_mu_prior(0,0,r,rr)) {
        sims(r,rr) = rnorm(trans_mu(0,0,r,rr), trans_mu_prior_sigma(0,0,r,rr));
        for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) sim_mu_prior_re(s,t,r,rr) = sims(r,rr);
      }
    }
    if((mu_model(r,rr) > 4) & (mu_model(r,rr) <= 8)) {//stock
      array<Type> sims(n_stocks, n_regions,n_regions-1);
      sims.setZero();
      for(int s = 0; s < n_stocks; s++){
        if(use_mu_prior(s,0,r,rr)) {
          sims(s,r,rr) = rnorm(trans_mu(s,0,r,rr), trans_mu_prior_sigma(s,0,r,rr));
          for(int t = 0; t < n_seasons; t++) sim_mu_prior_re(s,t,r,rr) = sims(s,r,rr);
        }
      }
    }
    if((mu_model(r,rr) > 8) & (mu_model(r,rr) <= 12)) {//season
      array<Type> sims(n_seasons,n_regions,n_regions-1);
      sims.setZero();
      for(int t = 0; t < n_seasons; t++){
        if(use_mu_prior(0,t,r,rr)) {
          sims(t,r,rr) = rnorm(trans_mu(0,t,r,rr), trans_mu_prior_sigma(0,t,r,rr));
          for(int s = 0; s < n_stocks; s++) sim_mu_prior_re(s,t,r,rr) = sims(t,r,rr);
        }
      }
    }
    if((mu_model(r,rr) > 12) & (mu_model(r,rr) <= 16)){ //stock,season
      for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++){
        if(use_mu_prior(s,t,r,rr)){
          sim_mu_prior_re(s,t,r,rr) = rnorm(trans_mu(s,t,r,rr), trans_mu_prior_sigma(s,t,r,rr));
        }
      }
    }
  }
  return(sim_mu_prior_re);
}
//done

template <class Type>
array<Type> get_nll_mu(array<Type> mu_repars, array<Type> mu_re, matrix<int> mu_model, array<int> can_move, vector<int> years_use){
  /* 
    get any nll components for time/age varying RE for movement parameters.
      mu_repars: n_stocks x n_seasons x n_regions x n_regions-1 x 3. parameters for distributions of random effects (sig, rho_a, rho_y)
          mu_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
                  mu_model: n_regions x n_regions-1. see definitions at top of move.hpp.
       can_move: n_stocks x n_seasons x n_regions x n_regions 0/1 whether fish can move from one region to another
      years_use: is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  //int n_y = mu_re.dim(3);
  int n_y = years_use.size();
  int n_regions = mu_re.dim(4);
  array<Type> nll(n_stocks,n_seasons,n_regions,n_regions-1);
  nll.setZero();

  array<int> can_move_reduced(n_stocks,n_seasons,n_regions,n_regions-1);
  for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
    int k = 0;
    for(int rr = 0; rr < n_regions; rr++) if(rr!=r) {
      can_move_reduced(s,t,r,k) = can_move(s,t,r,rr);
      k++;
    }
  }

  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++) {
    if((mu_model(r,rr) > 1) & (mu_model(r,rr) <=4)) if(can_move_reduced(0,0,r,rr)) {//constant, RE
      Type sigma_mu = exp(mu_repars(0,0,r,rr,0));
      Type rho_mu_a = geninvlogit(mu_repars(0,0,r,rr,1),Type(-1),Type(1),Type(1));
      Type rho_mu_y = geninvlogit(mu_repars(0,0,r,rr,2),Type(-1),Type(1),Type(1));
      if(mu_model(r,rr) == 2) { //age re
        Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
        vector<Type> mu_re_a(n_ages);
        for(int a = 0; a < n_ages; a++) mu_re_a(a) = mu_re(0,a,0,years_use(0),r,rr);
        nll(0,0,r,rr) += SCALE(AR1(rho_mu_a), Sigma_MU)(mu_re_a);
      }
      if(mu_model(r,rr) == 3) { //year re
        Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
        vector<Type> mu_re_y(n_y);
        for(int y = 0; y< n_y; y++) mu_re_y(y) = mu_re(0,0,0,years_use(y),r,rr);
        nll(0,0,r,rr) += SCALE(AR1(rho_mu_y), Sigma_MU)(mu_re_y);
      }
      if(mu_model(r,rr) == 4) { //age,year re
        Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
        array<Type> mu_re_ya(n_y,n_ages);
        for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) mu_re_ya(y,a) = mu_re(0,a,0,years_use(y),r,rr);
        nll(0,0,r,rr) += SCALE(SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)), Sigma_MU)(mu_re_ya); // must be array, not matrix!
      }
    }
    if((mu_model(r,rr) > 5) & (mu_model(r,rr) <=8)) {//stock, RE
      for(int s = 0; s < n_stocks; s++) if(can_move_reduced(s,0,r,rr)) {
        Type sigma_mu = exp(mu_repars(s,0,r,rr,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,0,r,rr,1),Type(-1),Type(1),Type(1));
        Type rho_mu_y = geninvlogit(mu_repars(s,0,r,rr,2),Type(-1),Type(1),Type(1));
        if(mu_model(r,rr) == 6) { //age re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
          vector<Type> mu_re_a(n_ages);
          for(int a = 0; a < n_ages; a++) mu_re_a(a) = mu_re(s,a,0,0,r,rr);
          nll(s,0,r,rr) += SCALE(AR1(rho_mu_a), Sigma_MU)(mu_re_a);
        }
        if(mu_model(r,rr) == 7) { //year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
          vector<Type> mu_re_y(n_y);
          for(int y = 0; y< n_y; y++) mu_re_y(y) = mu_re(s,0,0,years_use(y),r,rr);
          nll(s,0,r,rr) += SCALE(AR1(rho_mu_y), Sigma_MU)(mu_re_y);
        }
        if(mu_model(r,rr) == 8) { //age,year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
          array<Type> mu_re_ya(n_y,n_ages);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) mu_re_ya(y,a) = mu_re(s,a,0,years_use(y),r,rr);
          nll(s,0,r,rr) += SCALE(SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)), Sigma_MU)(mu_re_ya); // must be array, not matrix!
        }
      }
    }
    if((mu_model(r,rr) > 9) & (mu_model(r,rr) <=12)) {//season, RE
      for(int t = 0; t < n_seasons; t++) if(can_move_reduced(0,t,r,rr)){
        Type sigma_mu = exp(mu_repars(0,t,r,rr,0));
        Type rho_mu_a = geninvlogit(mu_repars(0,t,r,rr,1),Type(-1),Type(1),Type(1));
        Type rho_mu_y = geninvlogit(mu_repars(0,t,r,rr,2),Type(-1),Type(1),Type(1));
        if(mu_model(r,rr) == 10) { //age re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
          vector<Type> mu_re_a(n_ages);
          for(int a = 0; a < n_ages; a++) mu_re_a(a) = mu_re(0,a,t,0,r,rr);
          nll(0,t,r,rr) += SCALE(AR1(rho_mu_a), Sigma_MU)(mu_re_a);
        }
        if(mu_model(r,rr) == 11) { //year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
          vector<Type> mu_re_y(n_y);
          for(int y = 0; y< n_y; y++) mu_re_y(y) = mu_re(0,0,t,years_use(y),r,rr);
          nll(0,t,r,rr) += SCALE(AR1(rho_mu_y), Sigma_MU)(mu_re_y);
        }
        if(mu_model(r,rr) == 12) { //age,year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
          array<Type> mu_re_ya(n_y,n_ages);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) mu_re_ya(y,a) = mu_re(0,a,t,years_use(y),r,rr);
          nll(0,t,r,rr) += SCALE(SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)), Sigma_MU)(mu_re_ya); // must be array, not matrix!
        }
      }
    }
    if((mu_model(r,rr) > 13) & (mu_model(r,rr) <=16)) {//stock,season, RE
      for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) if(can_move_reduced(s,t,r,rr)) {
        Type sigma_mu = exp(mu_repars(s,t,r,rr,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,t,r,rr,1),Type(-1),Type(1),Type(1));
        Type rho_mu_y = geninvlogit(mu_repars(s,t,r,rr,2),Type(-1),Type(1),Type(1));
        if(mu_model(r,rr)  == 14) { //age re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
          vector<Type> mu_re_a(n_ages);
          for(int a = 0; a < n_ages; a++) mu_re_a(a) = mu_re(s,a,t,0,r,rr);
          nll(s,t,r,rr) += SCALE(AR1(rho_mu_a), Sigma_MU)(mu_re_a);
        }
        if(mu_model(r,rr)  == 15) { //year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
          vector<Type> mu_re_y(n_y);
          for(int y = 0; y< n_y; y++) mu_re_y(y) = mu_re(s,0,t,years_use(y),r,rr);
          nll(s,t,r,rr) += SCALE(AR1(rho_mu_y), Sigma_MU)(mu_re_y);
        }
        if(mu_model(r,rr)  == 16) { //age,year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
          array<Type> mu_re_ya(n_y,n_ages);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) mu_re_ya(y,a) = mu_re(s,a,t,years_use(y),r,rr);
          nll(s,t,r,rr) += SCALE(SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)), Sigma_MU)(mu_re_ya); // must be array, not matrix!
        }
      }
    }
  }
  return(nll);
}

template <class Type>
array<Type> simulate_mu_re(array<Type> mu_repars, array<Type> mu_re, matrix<int> mu_model, array<int> can_move, vector<int> years_use){
  /* 
    simulate andy time/age varying RE for movement parameters.
      mu_repars: n_stocks x n_seasons x n_regions x n_regions-1 x 3. parameters for distributions of random effects (sig, rho_a, rho_y)
          mu_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       mu_model: n_regions x n_regions-1. see definitions at top of move.hpp.
       can_move: n_stocks x n_seasons x n_regions x n_regions 0/1 whether fish can move from one region to another
      years_use: is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_y = years_use.size();
  //int n_y = mu_re.dim(3);
  int n_regions = mu_re.dim(4);
  array<Type> sim_mu_re = mu_re;//(n_stocks,n_ages,n_seasons,n_years_model, n_regions,n_regions-1);
  //sim_mu_re.setZero();

  array<int> can_move_reduced(n_stocks,n_seasons,n_regions,n_regions-1);
  for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
    int k = 0;
    for(int rr = 0; rr < n_regions; rr++) if(rr!=r) {
      can_move_reduced(s,t,r,k) = can_move(s,t,r,rr);
      k++;
    }
  }

  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++){
    if((mu_model(r,rr) > 1) & (mu_model(r,rr) <= 4)) if(can_move_reduced(0,0,r,rr)) {//constant, RE
      Type sigma_mu = exp(mu_repars(0,0,r,rr,0));
      Type rho_mu_a = geninvlogit(mu_repars(0,0,r,rr,1),Type(-1),Type(1),Type(1));
      Type rho_mu_y = geninvlogit(mu_repars(0,0,r,rr,2),Type(-1),Type(1),Type(1));
      if(mu_model(r,rr) == 2) { //age re
        Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
        vector<Type> mu_re_a(n_ages);
        AR1(rho_mu_a).simulate(mu_re_a);
        for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
          sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_a(a) * Sigma_MU;
        }
      }
      if(mu_model(r,rr) == 3) { //year re
        Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
        vector<Type> mu_re_y(n_y);
        AR1(rho_mu_y).simulate(mu_re_y);
        for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
          sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_y(y) * Sigma_MU;
        }
      }
      if(mu_model(r,rr) == 4) { //age,year re
        Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
        array<Type> mu_re_ya(n_y,n_ages);
        SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)).simulate(mu_re_ya);
        for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
          sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_ya(y,a) * Sigma_MU;
        }
      }
    }
    if((mu_model(r,rr) > 5) & (mu_model(r,rr) <= 8)) {//stock, RE
      for(int s = 0; s < n_stocks; s++) if(can_move_reduced(s,0,r,rr)){
        Type sigma_mu = exp(mu_repars(s,0,r,rr,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,0,r,rr,1),Type(-1),Type(1),Type(1));
        Type rho_mu_y = geninvlogit(mu_repars(s,0,r,rr,2),Type(-1),Type(1),Type(1));
        if(mu_model(r,rr) == 6) { //age re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
          vector<Type> mu_re_a(n_ages);
          AR1(rho_mu_a).simulate(mu_re_a);
          for(int t = 0; t < n_seasons; t++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
            sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_a(a) * Sigma_MU;
          }
        }
        if(mu_model(r,rr) == 7) { //year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
          vector<Type> mu_re_y(n_y);
          AR1(rho_mu_y).simulate(mu_re_y);
          for(int t = 0; t < n_seasons; t++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
            sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_y(y) * Sigma_MU;
          }
        }
        if(mu_model(r,rr) == 8) { //age,year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
          array<Type> mu_re_ya(n_y,n_ages);
          SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)).simulate(mu_re_ya);
          for(int t = 0; t < n_seasons; t++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
            sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_ya(y,a) * Sigma_MU;
          }
        }
      }
    }
    if((mu_model(r,rr) > 9) & (mu_model(r,rr) <= 12)) {//season, RE
      for(int t = 0; t < n_seasons; t++) if(can_move_reduced(0,t,r,rr)){ 
        Type sigma_mu = exp(mu_repars(0,t,r,rr,0));
        Type rho_mu_a = geninvlogit(mu_repars(0,t,r,rr,1),Type(-1),Type(1),Type(1));
        Type rho_mu_y = geninvlogit(mu_repars(0,t,r,rr,2),Type(-1),Type(1),Type(1));
        if(mu_model(r,rr) == 10) { //age re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
          vector<Type> mu_re_a(n_ages);
          AR1(rho_mu_a).simulate(mu_re_a);
          for(int s = 0; s < n_stocks; s++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) {
            sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_a(a) * Sigma_MU;
          }
        }
        if(mu_model(r,rr) == 11) { //year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
          vector<Type> mu_re_y(n_y);
          AR1(rho_mu_y).simulate(mu_re_y);
          for(int s = 0; s < n_stocks; s++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++){
            sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_y(y) * Sigma_MU;
          }
        }
        if(mu_model(r,rr) == 12) { //age,year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
          array<Type> mu_re_ya(n_y,n_ages);
          SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)).simulate(mu_re_ya);
          for(int s = 0; s < n_stocks; s++) for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++){
            sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_ya(y,a) * Sigma_MU;
          }
        }
      }
    }
    if((mu_model(r,rr) > 13) & (mu_model(r,rr) <= 16)) {//stock,season, RE
      for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) if(can_move_reduced(s,t,r,rr)){
        Type sigma_mu = exp(mu_repars(s,t,r,rr,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,t,r,rr,1),Type(-1),Type(1),Type(1));
        Type rho_mu_y = geninvlogit(mu_repars(s,t,r,rr,2),Type(-1),Type(1),Type(1));
        if(mu_model(r,rr) ==14) { //age re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_a,2)),-0.5); //marginal sd
          vector<Type> mu_re_a(n_ages);
          AR1(rho_mu_a).simulate(mu_re_a);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_a(a) * Sigma_MU;
        }
        if(mu_model(r,rr) ==15) { //year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)),-0.5); //marginal sd
          vector<Type> mu_re_y(n_y);
          AR1(rho_mu_y).simulate(mu_re_y);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++)  sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_y(y) * Sigma_MU;
        }
        if(mu_model(r,rr) ==16) { //age,year re
          Type Sigma_MU = sigma_mu * pow((1-pow(rho_mu_y,2)) * (1-pow(rho_mu_a,2)),-0.5); //marginal sd
          array<Type> mu_re_ya(n_y,n_ages);
          SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)).simulate(mu_re_ya);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++)  sim_mu_re(s,a,t,years_use(y),r,rr) = mu_re_ya(y,a) * Sigma_MU;
        }
      }
    }
  }
  return(sim_mu_re);
}
//done

template <class Type>
vector<Type> additive_ln_transform(vector<Type> x, int region, vector<int> can_move, int must_move){
  /* 
    use additive transformation (e.g., logistic-normal model)
    ensures that probabilities of moving and staying add to 1
              x: raw movement parameters 
         region: which region currently in
        do_move: which regions could move to
        can_move: 0/1 if can move from region to other regions 
      must_move: if 1, prob of leaving current region = 1

  */
  int D = x.size()+1;
  vector<Type> y(D);
  y.setZero();
  int j = 0;
  for(int i = 0; i < D; i++) {
    if(i != region) {
      if(can_move(i)) y(i) = exp(x(j)); //else prob of moving will be 0.
      j++;
    } 
    if(i == region) { //prob of staying will be 1- prob of moving
      if(must_move == 0) {
        y(i) = 1.0; //the last category/region (movement parameters not defined)
      } //else y(i) = 0.0 already. However, must make sure that trans_mu is fixed at 0 for one of the do_move regions to make returned y sum to 1.
    }
  }
  y /= sum(y);
  return(y);
}
//done

//provides transformed mu (good for sdreporting)
template<class Type>
array<Type> get_trans_mu_base(array<Type> trans_mu, array<Type>mu_re, array<Type> mu_prior_re, array<int> use_mu_prior, matrix<int> mu_model,
  array<Type> Ecov_lm, array<int> Ecov_how){
  /* 
    Construct base mu-at-age (excluding any density-dependence)
    currently continues random processes in any projection years!
          trans_mu: n_stocks x n_seasons x n_regions x n_regions-1 (mean) movement parameters
             mu_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       mu_prior_re: n_stocks x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
      use_mu_prior: n_stocks x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
          mu_model: n_regions x n_regions-1. see definitions at top of move.hpp.
           Ecov_lm: (n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_years_pop, n_Ecov) linear predictor for any Ecov effects on trans_mu_base
          Ecov_how: n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 values indicating to use effects on migration for each stock for each region (less 1).
  */
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_regions = mu_re.dim(4);
  int ny = mu_re.dim(3);
  //array<Type> Ecov_lm_mu(n_stocks, n_regions-1, n_ages, n_seasons, n_years_model + n_years_proj, n_Ecov);
  array<Type> trans_mu_base(n_stocks,n_ages,n_seasons,ny,n_regions, n_regions-1);
  trans_mu_base.setZero();
  
  if(n_regions>1){
    for(int r = 0; r< n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++) {
      for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y < ny; y++){
        if((mu_model(r,rr) > 0) & (mu_model(r,rr) <= 4)){ //constant
          if(mu_model(r,rr) == 2) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,a,0,0,r,rr); // age random effects
          if(mu_model(r,rr) == 3) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,0,0,y,r,rr); // year random effects
          if(mu_model(r,rr) == 4) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,a,0,y,r,rr); // age,year random effects
          if(use_mu_prior(0,0,r,rr)) trans_mu_base(s,a,t,y,r,rr) += mu_prior_re(0,0,r,rr);
          else trans_mu_base(s,a,t,y,r,rr) += trans_mu(s,t,r,rr);
        }
        if((mu_model(r,rr) > 4) & (mu_model(r,rr) <= 8)){ //stock
          if(mu_model(r,rr) == 6) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,a,0,0,r,rr); // age random effects
          if(mu_model(r,rr) == 7) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,0,0,y,r,rr); // year random effects
          if(mu_model(r,rr) == 8) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,a,0,y,r,rr); // age,year random effects
          if(use_mu_prior(s,0,r,rr)) trans_mu_base(s,a,t,y,r,rr) += mu_prior_re(s,0,r,rr);
          else trans_mu_base(s,a,t,y,r,rr) += trans_mu(s,t,r,rr);
        }
        if((mu_model(r,rr) > 8) & (mu_model(r,rr) <= 12)){ //season
          if(mu_model(r,rr) == 10) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,a,t,0,r,rr); // age random effects
          if(mu_model(r,rr) == 11) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,0,t,y,r,rr); // year random effects
          if(mu_model(r,rr) == 12) trans_mu_base(s,a,t,y,r,rr) += mu_re(0,a,t,y,r,rr); // age,year random effects
          if(use_mu_prior(0,t,r,rr)) trans_mu_base(s,a,t,y,r,rr) += mu_prior_re(0,t,r,rr);
        }
        if((mu_model(r,rr) > 12) & (mu_model(r,rr) <= 16)){ //stock,season
          if(mu_model(r,rr) == 14) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,a,t,0,r,rr); // age random effects
          if(mu_model(r,rr) == 15) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,0,t,y,r,rr); // year random effects
          if(mu_model(r,rr) == 16) trans_mu_base(s,a,t,y,r,rr) += mu_re(s,a,t,y,r,rr); // age,year random effects
          if(use_mu_prior(s,t,r,rr)) trans_mu_base(s,a,t,y,r,rr) += mu_prior_re(s,t,r,rr);
          else trans_mu_base(s,a,t,y,r,rr) += trans_mu(s,t,r,rr);
        }
        for(int i=0; i < Ecov_how.dim(0); i++) if(Ecov_how(i,s,a,t,r,rr) == 1) trans_mu_base(s,a,t,y,r,rr) += Ecov_lm(s,a,t,r,rr,y,i); //will be 0 if not used
      }
    }
  }
  //no projections options for mu. Just forecast any random or Ecov effects. Otherwise constant mu is the same as during model period.
  return(trans_mu_base); 
}
//done

template <class Type>
matrix<Type> get_mu_matrix(int stock, int age, int season, int year, vector<int> mig_type, array<int> can_move, array<int> must_move, array<Type> trans_mu_base){
  /* 
    Construct n_regions x n_regions movement matrix
      stock: which stock
      age: which age
      season: which season
      year: which year
      mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
      can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
      must_move: n_stocks x n_seasons x n_regions: 0/1 determining if it must leave the region
      trans_mu_base: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions-1. array retruned by get_trans_mu_base
  */

  int n_regions = trans_mu_base.dim(4);
  matrix<Type> mu(n_regions,n_regions);
  mu.setZero();
  if(mig_type(stock) == 0) //migration is instantaneous after survival and mortality, so P is easy.
  {
    for(int r = 0; r < n_regions; r++) {
      vector<Type> trans_par(n_regions-1);
      vector<int> can_move_r(n_regions); //can_move_r(r) is ignored because dictated by must_move?
      for(int rr = 0; rr < n_regions-1; rr++) trans_par(rr) = trans_mu_base(stock,age,season,year,r,rr);
      for(int j = 0; j < n_regions; j++) can_move_r(j) = can_move(stock,season,r,j);
      vector<Type> pmove = additive_ln_transform(trans_par, r, can_move_r, must_move(stock,season,r));
      for(int j = 0; j < n_regions; j++) mu(r,j) = pmove(j);
    }
  }
  if(mig_type(stock) == 1) //migration occurs continuously during interval, return infinitesimal generator.
  {
    for(int i = 0; i < n_regions; i++) {
      int k = 0;
      for(int j = 0; j < n_regions; j++){ 
        if(j!=i) {
          if(can_move(stock,season,i,j)) mu(i,j) = exp(trans_mu_base(stock,age,season,year,i,k)); //log of transition intensities
          k++; //max k = n_regions -1 (-1)
        }
      }
    }
    for(int r = 0; r< n_regions; r++) mu(r,r) = -(mu.row(r)).sum(); //hazard
  }
  return(mu);
}
//done

template <class Type>
array<Type> get_avg_mu(array<Type> trans_mu_base, vector<int> years, vector<int> mig_type, array<int> can_move,
 array<int> must_move){
  /* 
    Construct n_stocks x n_ages x n_seasons x n_regions x n_regions array of "averaged" movement parameters over years
      stock: which stock
      age: which age
      season: which season
      years: which years to average over
      mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
      can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
      must_move: n_stocks x n_seasons x n_regions: 0/1 determining if it must leave the region
      trans_mu_base: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions-1. array retruned by get_trans_mu_base
  */

  int n_stocks =  trans_mu_base.dim(0);
  int n_ages = trans_mu_base.dim(1);
  int n_y = years.size();
  int n_seasons = trans_mu_base.dim(2);
  int n_regions = trans_mu_base.dim(4);
  array<Type> avg_mu(n_stocks, n_ages,n_seasons,n_regions,n_regions);
  avg_mu.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a ++) for(int t = 0; t < n_seasons; t++){
    //matrix<Type> mu(n_regions,n_regions);
    //mu.setZero();
    if(mig_type(s) == 0) { //migration is instantaneous after survival and mortality, so P is easy.
      // from each region, average the probabilities of movement and then rescale. 
      for(int r = 0; r < n_regions; r++) {
        vector<Type> trans_par(n_regions-1);
        trans_par.setZero();
        vector<int> can_move_r(n_regions); //can_move_r(r) is ignored because dictated by must_move?
        can_move_r.setZero();
        vector<Type> pmove(n_regions);
        pmove.setZero();
        for(int j = 0; j < n_regions; j++) can_move_r(j) = can_move(s,t,r,j);
        for(int y = 0; y < n_y; y++) {
          for(int rr = 0; rr < n_regions-1; rr++) trans_par(rr) = trans_mu_base(s,a,t,years(y),r,rr);
          pmove += additive_ln_transform(trans_par, r, can_move_r, must_move(s,t,r))/Type(n_y);
        }
        //pmove /= sum(pmove); //NOT needed to ensure averaged values sum to 1
        for(int j = 0; j < n_regions; j++) avg_mu(s,a,t,r,j) = pmove(j);
      }
    }
    if(mig_type(s) == 1) { //migration occurs continuously during interval, return infinitesimal generator.
      //for each region, average the yearly instantaneous movement rates to other regions (analogous to M and F)
      for(int r = 0; r < n_regions; r++) {
        int k = 0;
        for(int j = 0; j < n_regions; j++){ 
          if(j!=r) {
            if(can_move(s,t,r,j)) for(int y = 0; y < n_y; y++) {
              avg_mu(s,a,t,r,j) = exp(trans_mu_base(s,a,t,years(y),r,k))/Type(n_y); //log of transition intensities
              
            }
            avg_mu(s,a,t,r,r) -= avg_mu(s,a,t,r,j); //-sum = hazard
            k++; //max k = n_regions -1 (-1)
          }
        }
      }
    }
  }
  return avg_mu;
}

//all movement matrices
template <class Type>
array<Type> get_mu(array<Type> trans_mu_base, array<int> can_move,  array<int> must_move, vector<int> mig_type, 
  int n_years_proj, int n_years_model, int proj_mu_opt, vector<int> avg_years){
  /* 
    Construct n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
      trans_mu_base: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions-1. array retruned by get_trans_mu_base
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
          must_move: n_stocks x n_seasons x n_regions: 0/1 determining if it must leave the region
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
       n_years_proj: number of projection years
      n_years_model: number of years before projections
        proj_mu_opt: 1: use averega of mu over years, 2: use random trans_mu_base with RE and/or Ecov effects in projection years
          avg_years: which model years to use for averaging mu
  */
  int n_stocks = trans_mu_base.dim(0);
  int n_ages = trans_mu_base.dim(1);
  int n_seasons = trans_mu_base.dim(2);
  int n_y = trans_mu_base.dim(3);
  int n_regions = trans_mu_base.dim(4);
  array<Type> mu(n_stocks,n_ages, n_seasons, n_y, n_regions,n_regions);
  mu.setZero();
  if(n_regions>1) {
    for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y < n_y; y++){
      matrix<Type> mu_y = get_mu_matrix(s,a,t,y,mig_type,can_move,must_move,trans_mu_base);
      for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) mu(s,a,t,y,r,rr) = mu_y(r,rr);
    }
    // add to mu in projection years
    if(n_years_proj>0){ 
      //int n_toavg = avg_years.size();
      if(proj_mu_opt == 2){ // use average mu over avg.yrs 
        array<Type> mu_avg = get_avg_mu(trans_mu_base,avg_years,mig_type, can_move, must_move);
        for(int y = n_years_model; y < n_y; y++)for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) {
          for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
            mu(s,a,t,y,r,rr) = mu_avg(s,a,t,r,rr);
          }
        }
      } else { // proj_mu_opt == 1, use mu_re and/or ecov_lm in projection years
        for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++){ 
          for(int y = n_years_model; y < n_y; y++){
            matrix<Type> mu_y = get_mu_matrix(s,a,t,y,mig_type,can_move,must_move,trans_mu_base);
            for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) mu(s,a,t,y,r,rr) = mu_y(r,rr);
          }
        }
      }
    }
  }
  return(mu);
}
//done

//extract array of mu parameters for a given year
template <class Type>
array<Type> get_mu_y(int y, array<Type> mu){
  int n_stocks = mu.dim(0);
  int n_ages = mu.dim(1);
  int n_seasons = mu.dim(2);
  //int n_y = trans_mu_base.dim(3);
  int n_regions = mu.dim(4);
  array<Type> mu_y(n_stocks,n_ages, n_seasons, n_regions,n_regions);
  mu_y.setZero();

  if(n_regions>1) {
    for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++){
      for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++)  mu_y(s,a,t,r,rr) = mu(s,a,t,y,r,rr);
    }
  }
  return mu_y;
}

//sdreport a reduced set of trans_mu_base parameters based on mu_model
template <class Type>
array<int> get_mu_sdrep_indices(matrix<int> mu_model, array<Type> trans_mu_base){
  int n_stocks = trans_mu_base.dim(0);
  int n_ages = trans_mu_base.dim(1);
  int n_seasons = trans_mu_base.dim(2);
  int n_y = trans_mu_base.dim(3);
  int n_regions = trans_mu_base.dim(4);
  array<int> mu_index(n_regions,n_regions-1,2);
  mu_index.setZero();
  int k = 0;
  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++) {
    mu_index(r,rr,0) = k;
    if(mu_model(r,rr) == 1) k++; //do nothing 0,0 already, increment k
    if(mu_model(r,rr) == 2) { 
      mu_index(r,rr,1) = k + n_ages-1;
      k += n_ages;
    }
    if(mu_model(r,rr) == 3) { 
      mu_index(r,rr,1) = k + n_y-1;
      k += n_y;
    }
    if(mu_model(r,rr) == 4) { 
      mu_index(r,rr,1) = k + n_ages * n_y - 1;
      k += n_ages * n_y;
    }
    if(mu_model(r,rr) == 5) { 
      mu_index(r,rr,1) = k + n_stocks-1;
      k += n_stocks;
    }
    if(mu_model(r,rr) == 6) { 
      mu_index(r,rr,1) = k + n_stocks*n_ages-1;
      k += n_stocks*n_ages;
    }
    if(mu_model(r,rr) == 7) { 
      mu_index(r,rr,1) = k + n_stocks*n_y-1;
      k += n_stocks*n_y;
    }
    if(mu_model(r,rr) == 8) { 
      mu_index(r,rr,1) = k + n_stocks*n_y*n_ages-1;
      k += n_stocks*n_y*n_ages;
    }
    if(mu_model(r,rr) == 9) { 
      mu_index(r,rr,1) = k + n_seasons-1;
      k += n_seasons;
    }
    if(mu_model(r,rr) == 10) { 
      mu_index(r,rr,1) = k + n_seasons*n_ages-1;
      k += n_seasons*n_ages;
    }
    if(mu_model(r,rr) == 11) { 
      mu_index(r,rr,1) = k + n_seasons*n_y-1;
      k += n_seasons*n_y;
    }
    if(mu_model(r,rr) == 12) { 
      mu_index(r,rr,1) = k + n_seasons*n_y*n_ages-1;
      k += n_seasons*n_y*n_ages;
    }
    if(mu_model(r,rr) == 13) { 
      mu_index(r,rr,1) = k + n_stocks*n_seasons-1;
      k += n_stocks*n_seasons;
    }
    if(mu_model(r,rr) == 14) { 
      mu_index(r,rr,1) = k + n_stocks*n_seasons*n_ages-1;
      k += n_stocks*n_seasons*n_ages;
    }
    if(mu_model(r,rr) == 15) { 
      mu_index(r,rr,1) = k + n_stocks*n_seasons*n_y-1;
      k += n_stocks*n_seasons*n_y;
    }
    if(mu_model(r,rr) == 16) { 
      mu_index(r,rr,1) = k + n_stocks*n_seasons*n_y*n_ages-1;
      k += n_stocks*n_seasons*n_y*n_ages;
    }
  }
  return mu_index;
}

template <class Type>
vector <Type> get_trans_mu_base_sdrep(array<Type> trans_mu_base, matrix<int> mu_model, array<int> mu_sdrep_index){
  int n_stocks = trans_mu_base.dim(0);
  int n_ages = trans_mu_base.dim(1);
  int n_seasons = trans_mu_base.dim(2);
  int n_y = trans_mu_base.dim(3);
  int n_regions = trans_mu_base.dim(4);
  int sdrep_size = mu_sdrep_index(n_regions-1,n_regions-2,1)+1; //dim = n_regions x n_regions -1 x 2
  vector<Type> trans_mu_sdrep(sdrep_size);
  trans_mu_sdrep.setZero();
  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions-1; rr++) {
    int ind = mu_sdrep_index(r,rr,0);
    if(mu_model(r,rr) == 1){ //constant
      trans_mu_sdrep(ind) = trans_mu_base(0,0,0,0,r,rr);
    }
    if(mu_model(r,rr) == 2){ //age random effects
      for(int a = 0; a < n_ages; a++){  
        trans_mu_sdrep(ind+a) = trans_mu_base(0,a,0,0,r,rr);
      }
    }
    if(mu_model(r,rr) == 3){ //year random effects
      for(int y = 0; y < n_y; y++) {
        trans_mu_sdrep(ind+y) = trans_mu_base(0,0,0,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 4){ //age,year random effects
      for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_y; y++){ 
        trans_mu_sdrep(ind + y + n_y * a) = trans_mu_base(0,a,0,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 5){ //stock
      for(int s = 0; s< n_stocks; s++) {
        trans_mu_sdrep(ind + s) = trans_mu_base(s,0,0,0,r,rr);
      }
    }
    if(mu_model(r,rr) == 6){ //stock, age random effects
      for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) {
        trans_mu_sdrep(ind + a + n_ages * s) = trans_mu_base(s,a,0,0,r,rr);
      }
    }
    if(mu_model(r,rr) == 7){ //stock, year random effects
      for(int s = 0; s< n_stocks; s++) for(int y = 0; y < n_y; y++) {
        trans_mu_sdrep(ind + y + n_y * s) = trans_mu_base(s,0,0,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 8){ //stock, age,year random effects
      for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_y; y++) {
        trans_mu_sdrep(ind + y + n_y * (a + n_ages * s)) = trans_mu_base(s,a,0,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 9){ //season
      for(int t = 0; t < n_seasons; t++) trans_mu_sdrep(ind + t) = trans_mu_base(0,0,t,0,r,rr);      
    }
    if(mu_model(r,rr) == 10){ //season, age random effects
      for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) {
        trans_mu_sdrep(ind + a + n_ages * t) = trans_mu_base(0,a,t,0,r,rr);
      }
    }
    if(mu_model(r,rr) == 11){ //season, year random effects
      for(int t = 0; t < n_seasons; t++) for(int y = 0; y < n_y; y++) {
        trans_mu_sdrep(ind + y + n_y * t) = trans_mu_base(0,0,t,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 12){ //season, age,year random effects
      for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) for(int y = 0; y < n_y; y++) {
        trans_mu_sdrep(ind + y + n_y * (a + n_ages * t)) = trans_mu_base(0,a,t,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 13){ //stock,season
      for(int s = 0; s< n_stocks; s++) for(int t = 0; t < n_seasons; t++) {
        trans_mu_sdrep(ind + t + n_seasons*t) = trans_mu_base(s,0,t,0,r,rr);
      }
    }
    if(mu_model(r,rr) == 14){ //stock,season, age random effects
      for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) {
        trans_mu_sdrep(ind + a + n_ages * (t + n_seasons * s)) = trans_mu_base(s,a,t,0,r,rr);
      }
    }
    if(mu_model(r,rr) == 15){ //stock,season, year random effects
      for(int s = 0; s< n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y < n_y; y++) {
        trans_mu_sdrep(ind + y + n_y * (t + n_seasons*s)) = trans_mu_base(s,0,t,y,r,rr);
      }
    }
    if(mu_model(r,rr) == 16){ //stock,season, year random effects
      for(int s = 0; s< n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) {
        trans_mu_sdrep(ind + y + n_y * (a + n_ages * (t + n_seasons * s))) = trans_mu_base(s,a,t,y,r,rr);
      }
    }
  }
  return trans_mu_sdrep;
}
