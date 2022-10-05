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
array<Type> get_nll_mu_prior(array<Type> mu_prior_re, array<Type> trans_mu, array<Type> trans_mu_prior_sigma, array<int> use_mu_prior, int mu_model){
  /* 
    get any nll components for priors/posteriors on mu parameters.
               mu_prior_re: n_stocks x n_ages x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
                  trans_mu: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
      trans_mu_prior_sigma: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
              use_mu_prior: n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
                  mu_model: see definitions at top of move.hpp.
  */
  int n_stocks = trans_mu.dim(0);
  int n_ages = trans_mu.dim(1);
  int n_seasons = trans_mu.dim(2);
  int n_regions = trans_mu.dim(3);
  array<Type> nll(n_stocks,n_ages, n_seasons, n_regions, n_regions-1);
  nll.setZero();
  if((mu_model == 1) | (mu_model == 2) | (mu_model == 3) | (mu_model == 4) ){
    for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(0,0,t,r,rr)) nll(0,0,t,r,rr) -= dnorm(mu_prior_re(0,0,t,r,rr), trans_mu(0,0,t,r,rr), trans_mu_prior_sigma(0,0,t,r,rr),1);
    }
  }
  if((mu_model == 5) | (mu_model == 6) | (mu_model == 7) | (mu_model == 8)){
    for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(s,0,t,r,rr)) nll(s,0,t,r,rr) -= dnorm(mu_prior_re(s,0,t,r,rr), trans_mu(s,0,t,r,rr), trans_mu_prior_sigma(s,0,t,r,rr),1);
    }
  }
  return(nll);
}
//done

template <class Type>
array<Type> simulate_mu_prior_re(array<Type> mu_prior_re, array<Type> trans_mu, array<Type> trans_mu_prior_sigma, array<int> use_mu_prior, int mu_model){
  /* 
    simulate and RE for priors on mu parameters.
               mu_prior_re: n_stocks x n_ages x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
                  trans_mu: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
      trans_mu_prior_sigma: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
              use_mu_prior: n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
                  mu_model: see definitions at top of move.hpp.
  */
  int n_stocks = trans_mu.dim(0);
  int n_ages = trans_mu.dim(1);
  int n_seasons = trans_mu.dim(2);
  int n_regions = trans_mu.dim(3);
  array<Type> sim_mu_prior_re(n_stocks,n_ages, n_seasons, n_regions, n_regions-1);
  sim_mu_prior_re.setZero();
  if((mu_model == 1) | (mu_model == 2) | (mu_model == 3) | (mu_model == 4) ){
    for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(0,0,t,r,rr)==1) {
        sim_mu_prior_re(0,0,t,r,rr) = rnorm(trans_mu(0,0,t,r,rr), trans_mu_prior_sigma(0,0,t,r,rr));
      }
    }
  }
  if((mu_model == 5) | (mu_model == 6) | (mu_model == 7) | (mu_model == 8)){
    for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr< n_regions-1; rr++){
      if(use_mu_prior(s,0,t,r,rr)){
        sim_mu_prior_re(s,0,t,r,rr) = rnorm(trans_mu(s,0,t,r,rr), trans_mu_prior_sigma(s,0,t,r,rr));
      }
    }
  }
  return(sim_mu_prior_re);
}
//done

template <class Type>
array<Type> get_nll_mu(array<Type> mu_repars, array<Type> mu_re, int mu_model){
  /* 
    get any nll components for time/age varying RE for movement parameters.
      mu_repars: n_stocks x n_ages x n_seasons x n_regions x 4. parameters for distributions of random effects (sig, rho_a, rho_y, rho_r)
          mu_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       mu_model: see definitions at top of move.hpp.
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_y = mu_re.dim(3);
  int n_regions = mu_re.dim(4);
  array<Type> nll(n_stocks,n_ages,n_seasons,n_regions,n_regions-1);
  nll.setZero();
  if((mu_model != 1) & (mu_model != 5)){ //some type of random effects
    int ns = 1;
    if(mu_model > 5) ns = n_stocks; 
    if((mu_model == 2) | (mu_model == 6)){ // age random effects
      for(int s = 0; s < ns; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++){
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,0,t,r,1),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = sigma_mu / pow((1-pow(rho_mu_a,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_a(n_ages);
          for(int a = 0; a < n_ages; a++) mu_re_a(a) = mu_re(s,a,t,0,r,rr);
          nll(s,0,t,r,rr) += SCALE(AR1(rho_mu_a), Sigma_MU)(mu_re_a);
        }
      }
    }
    if((mu_model == 3) | (mu_model == 7)){ // year random effects
      for(int s = 0; s < ns; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_y = geninvlogit(mu_repars(s,0,t,r,2),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = sigma_mu / pow((1-pow(rho_mu_y,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_y(n_y);
          for(int y = 0; y< n_y; y++) mu_re_y(y) = mu_re(s,0,t,y,r,rr);
          nll(s,0,t,r,rr) += SCALE(AR1(rho_mu_y), Sigma_MU)(mu_re_y);
        }
      }
    }
    if((mu_model == 4) | (mu_model == 8)){ // age and year random effects
      for(int s = 0; s < ns; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,0,t,r,1),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type rho_mu_y = geninvlogit(mu_repars(s,0,t,r,2),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = sigma_mu * pow(((1-pow(rho_mu_y,2))*(1-pow(rho_mu_a,2))),-0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          array<Type> mu_re_ya(n_y,n_ages);
          for(int y = 0; y< n_y; y++) for(int a = 0; a < n_ages; a++) mu_re_ya(y,a) = mu_re(s,a,t,y,r,rr);
          nll(s,0,t,r,rr) += SCALE(SEPARABLE(AR1(rho_mu_a),AR1(rho_mu_y)), Sigma_MU)(mu_re_ya); // must be array, not matrix!
        }
      }
    }
  }
  return(nll);
}
//done

template <class Type>
array<Type> simulate_mu_re(array<Type> mu_repars, array<Type> mu_re, int mu_model){
  /* 
    simulate andy time/age varying RE for movement parameters.
      mu_repars: n_stocks x n_ages x n_seasons x n_regions x 4. parameters for distributions of random effects (sig, rho_a, rho_y, rho_r)
          mu_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       mu_model: see definitions at top of move.hpp.
  */
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_y = mu_re.dim(3);
  int n_regions = mu_re.dim(4);
  array<Type> sim_mu_re(n_stocks,n_ages,n_seasons,n_regions,n_regions-1);
  sim_mu_re.setZero();
  if((mu_model != 1) & (mu_model != 5)){ //some type of random effects
    int ns = 1;
    if(mu_model > 5) ns = n_stocks; 
    if((mu_model == 2) | (mu_model == 6)){ // age random effects
      for(int s = 0; s < ns; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++){
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,0,t,r,1),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_a,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_a(n_ages);
          AR1(rho_mu_a).simulate(mu_re_a);
          for(int a = 0; a < n_ages; a++) sim_mu_re(s,a,t,0,r,rr) = mu_re_a(a) * Sigma_MU;
        }
      }
    }
    if((mu_model == 3) | (mu_model == 7)){ // year random effects
      for(int s = 0; s < ns; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_y = geninvlogit(mu_repars(s,0,t,r,2),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type Sigma_MU = pow(pow(sigma_mu,2) / (1-pow(rho_mu_y,2)),0.5); //marginal sd
        for(int rr = 0; rr < n_regions-1; rr++){
          vector<Type> mu_re_y(n_y);
          AR1(rho_mu_y).simulate(mu_re_y);
          for(int y = 0; y< n_y; y++) sim_mu_re(s,0,t,y,r,rr) = mu_re_y(y) * Sigma_MU;
        }
      }
    }
    if((mu_model == 4) | (mu_model == 8)){ // age and year random effects
      for(int s = 0; s < ns; s++) for(int t = 0; t < n_seasons; t++) for(int r = 0; r < n_regions; r++) {
        Type sigma_mu = exp(mu_repars(s,0,t,r,0));
        Type rho_mu_a = geninvlogit(mu_repars(s,0,t,r,1),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
        Type rho_mu_y = geninvlogit(mu_repars(s,0,t,r,2),Type(-1),Type(1),Type(1));//rho_trans(mu_repars(s,0,t,r,2));
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
//done

template <class Type>
vector<Type> additive_ln_transform(vector<Type> x, int region, vector<int> can_move, int must_move){
  /* 
    use additive transformation (e.g., logistic-normal model)
    ensures that probabilities of moving and staying add to 1
              x: raw movement parameters 
         region: which region currently in
        do_move: which regions could move to
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
array<Type> get_trans_mu_base(array<Type> trans_mu, array<Type>mu_re, array<Type> mu_prior_re, array<int> use_mu_prior, int mu_model,
  array<Type> Ecov_lm, array<int> use_Ecov){
  /* 
    Construct base mu-at-age (excluding any density-dependence)
    currently continues random processes in any projection years!
          trans_mu: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 (mean) movement parameters
             mu_re: n_stocks x n_ages x n_seasons x n_y x n_regions x n_regions-1. RE for movement.
       mu_prior_re: n_stocks x n_ages x n_seasons x n_regions x n_regions-1. RE for posterior (given prior) (mean) movement parameters
      use_mu_prior: n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 whether to apply prior for each movement parameter
          mu_model: see definitions at top of move.hpp.
           Ecov_lm: (n_stocks, n_ages, n_seasons, n_regions, n_regions-1, n_years_pop, n_Ecov) linear predictor for any Ecov effects on trans_mu_base
          use_Ecov: n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions-1: 0/1 values indicating to use effects on migration for each stock for each region (less 1).
  */
  int n_stocks = mu_re.dim(0);
  int n_ages = mu_re.dim(1);
  int n_seasons = mu_re.dim(2);
  int n_regions = mu_re.dim(4);
  int ny = mu_re.dim(3);
  //array<Type> Ecov_lm_mu(n_stocks, n_regions-1, n_ages, n_seasons, n_years_model + n_years_proj, n_Ecov);
  array<Type> trans_mu_base(n_stocks,n_ages,n_seasons,ny,n_regions, n_regions-1);
  trans_mu_base.setZero();
  
  if(n_regions>1) for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y < ny; y++){
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
      for(int i=0; i < use_Ecov.dim(0); i++) if(use_Ecov(i,s,a,t,r,rr) == 1) trans_mu_base(s,a,t,y,r,rr) += Ecov_lm(s,a,t,r,rr,y,i); //will be 0 if not used
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
      can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
      must_move: n_stocks x ages x n_seasons x n_regions: 0/1 determining if it must leave the region
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
      for(int j = 0; j < n_regions; j++) can_move_r(j) = can_move(stock,age,season,r,j);
      vector<Type> pmove = additive_ln_transform(trans_par, r, can_move_r, must_move(stock,age,season,r));
      for(int j = 0; j < n_regions; j++) mu(r,j) = pmove(j);
    }
  }
  if(mig_type(stock) == 1) //migration occurs continuously during interval, return infinitesimal generator.
  {
    for(int i = 0; i < n_regions; i++) {
      int k = 0;
      for(int j = 0; j < n_regions; j++){ 
        if(j!=i) {
          k++; //max k = n_regions -1 (-1)
          if(can_move(i,j)==1) mu(i,j) = exp(trans_mu_base(stock,age,season,year,i,k)); //log of transition intensities
        }
      }
    }
    for(int r = 0; r< n_regions; r++) mu(r,r) = -(mu.row(r)).sum(); //hazard
  }
  return(mu);
}
//done

//all movement matrices
template <class Type>
array<Type> get_mu(array<Type> trans_mu_base, array<int> can_move,  array<int> must_move, vector<int> mig_type){
  /* 
    Construct n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
      trans_mu_base: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions-1. array retruned by get_trans_mu_base
      can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
      must_move: n_stocks x ages x n_seasons x n_regions: 0/1 determining if it must leave the region
      mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
  */
  int n_stocks = trans_mu_base.dim(0);
  int n_ages = trans_mu_base.dim(1);
  int n_seasons = trans_mu_base.dim(2);
  int n_y = trans_mu_base.dim(3);
  int n_regions = trans_mu_base.dim(4);
  array<Type> mu(n_stocks,n_ages, n_seasons, n_y, n_regions,n_regions);
  mu.setZero();
  if(n_regions>1) for(int s = 0; s< n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int t = 0; t < n_seasons; t++) for(int y = 0; y < n_y; y++){
    matrix<Type> mu_y = get_mu_matrix(s,a,t,y,mig_type,can_move,must_move,trans_mu_base);
    for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) mu(s,a,t,y,r,rr) = mu_y(r,rr);
  }
  return(mu);
}
//done


template <class Type>
matrix<Type> get_avg_mu_matrix(int stock, int age, int season, vector<int> years, vector<int> mig_type, array<int> can_move,
 array<int> must_move, array<Type> trans_mu_base){
  /* 
    Construct n_regions x n_regions movement matrix "averaged" over years
      stock: which stock
      age: which age
      season: which season
      years: which years to average over
      mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
      can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
      must_move: n_stocks x ages x n_seasons x n_regions: 0/1 determining if it must leave the region
      trans_mu_base: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions-1. array retruned by get_trans_mu_base
  */

  int n_regions = trans_mu_base.dim(4);
  matrix<Type> mu(n_regions,n_regions);
  mu.setZero();
  if(mig_type(stock) == 0) { //migration is instantaneous after survival and mortality, so P is easy.
    // from each region, average the probabilities of movement and then rescale. 
    for(int r = 0; r < n_regions; r++) {
      vector<Type> trans_par(n_regions-1);
      trans_par.setZero();
      vector<int> can_move_r(n_regions); //can_move_r(r) is ignored because dictated by must_move?
      can_move_r.setZero();
      vector<Type> pmove(n_regions);
      pmove.setZero();
      for(int j = 0; j < n_regions; j++) can_move_r(j) = can_move(stock,age,season,r,j);
      for(int y = 0; y < years.size(); y++) {
        for(int rr = 0; rr < n_regions-1; rr++) trans_par(rr) = trans_mu_base(stock,age,season,years(y),r,rr);
        pmove += additive_ln_transform(trans_par, r, can_move_r, must_move(stock,age,season,r))/years.size();
      }
      pmove = pmove/sum(pmove);
      for(int j = 0; j < n_regions; j++) mu(r,j) = pmove(j);
    }
  }
  if(mig_type(stock) == 1) { //migration occurs continuously during interval, return infinitesimal generator.
    //for each region, average the yearly instantaneous movement rates to other regions (analogous to M and F)
    for(int i = 0; i < n_regions; i++) {
      int k = 0;
      for(int j = 0; j < n_regions; j++){ 
        if(j!=i) {
          k++; //max k = n_regions -1 (-1)
          if(can_move(i,j)==1) for(int y = 0; y < years.size(); y++) {
            mu(i,j) = exp(trans_mu_base(stock,age,season,years(y),i,k))/years.size(); //log of transition intensities
          }
        }
      }
    }
    for(int r = 0; r< n_regions; r++) mu(r,r) = -(mu.row(r)).sum(); //hazard
  }
  return(mu);
}
//done
