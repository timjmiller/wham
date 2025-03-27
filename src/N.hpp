
template <class T>
array<T> get_marginal_NAA_sigma(array<T> log_NAA_sigma, array<T> trans_NAA_rho, vector<int> NAA_re_model, int decouple_recruitment = 0){

  int n_stocks = log_NAA_sigma.dim(0);
  int n_regions = log_NAA_sigma.dim(1);
  int n_ages = log_NAA_sigma.dim(2);
  int n_rho = trans_NAA_rho.dim(2);
  array<T> NAA_rho = trans_NAA_rho;
  array<T> NAA_sigma = log_NAA_sigma;
  array<T> marginal_sigma = log_NAA_sigma;
  marginal_sigma.setZero();
  int rho_y_R_ind = 1;
  if(decouple_recruitment) rho_y_R_ind = 2;
  T NAA_rho_y = 0, NAA_rho_a = 0;;

  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
    for(int k = 0; k < n_rho ; k++) NAA_rho(s,r,k) = geninvlogit(trans_NAA_rho(s,r,k), T(-1), T(1), T(1)); //using scale =1 ,2 is legacy
    if((NAA_re_model(s) == 1) | ((NAA_re_model(s) == 2) & decouple_recruitment)){ //"rec"
      NAA_rho_y = NAA_rho(s,r,rho_y_R_ind);
      marginal_sigma(s,r,0) = NAA_sigma(s,r,0) * pow(1-pow(NAA_rho_y,2),-0.5);
    }
    if(NAA_re_model(s) == 2){
      NAA_rho_y = geninvlogit(trans_NAA_rho(s,r,1), T(-1), T(1), T(1)); //using scale =1 ,2 is legacy        
      NAA_rho_a = geninvlogit(trans_NAA_rho(s,r,0), T(-1), T(1), T(1)); //using scale =1 ,2 is legacy
      int age_start = 0;
      if(decouple_recruitment) age_start = 1;
      for(int a = age_start; a < n_ages ; a++) {
        NAA_sigma(s,r,a) = exp(log_NAA_sigma(s,r,a));
        marginal_sigma(s,r,a) = NAA_sigma(s,r,a) * pow((1-pow(NAA_rho_y,2))*(1-pow(NAA_rho_a,2)),-0.5);
      }
    }
  }
  return(marginal_sigma);
}


template <class Type>
matrix<Type> get_nll_N1(vector<int> N1_model, array<Type>log_N1, array<Type> N1_repars, array<int> NAA_where) {
  /* 
    get nll contribtions for any N1 random effects
       N1_model: 0: (n_stocks) just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
         log_N1: (n_stocks x n_regions x n_ages) fixed or random effects for initial numbers at age
      N1_repars: (n_stocks x 3) mean, sig, rho; sd and correlation parameters for N1 random effects
      NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
  */
  //AR1 RE for N1 if N1_model = 2, mapped appropriately on R side
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(2);
  matrix<Type> nll(n_stocks,n_regions);
  nll.setZero();
  for(int s = 0; s < n_stocks; s++) if(N1_model(s) ==2) for(int r = 0; r < n_regions; r++) {
    int n_ages_r = 0;
    //need to count how many age classes exist in region r on Jan 1.
    for(int a = 0; a < n_ages; a++) if(NAA_where(s,r,a)) n_ages_r++; 
    if(n_ages_r>0){
      vector<Type> re_sr(n_ages_r);
      Type mu = N1_repars(s,r,0);
      Type rho = geninvlogit(N1_repars(s,r,2),Type(-1),Type(1),Type(1));
      Type sigma = exp(N1_repars(s,r,1)) * pow(1 - pow(rho,2),-0.5); //marginal variance
      re_sr.setZero();
      int a_count = 0;
      for(int a = 0; a < n_ages; a++) {
        if(NAA_where(s,r,a)) {
          re_sr(a_count) = log_N1(s,r,a) - mu;
          a_count++;
        } 
      }
      nll(s,r) += SCALE(AR1(rho), sigma)(re_sr);
    }
  }
  return(nll);
}
//Done.

template <class Type>
array<Type> simulate_log_N1(vector<int> N1_model, array<Type>log_N1, array<Type> N1_repars, array<int> NAA_where)
{ 
  /* 
    simulate any N1 random effects
       N1_model: 0: (n_stocks) just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
         log_N1: (n_stocks x n_ages) current array of random effects for initial numbers at age (used for size information)
      N1_repars: (n_stocks x 3) mean, sig, rho; sd and correlation parameters for N1 random effects
      NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
  */
  //only used if N1_model = 2
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(2);
  array<Type> sim_log_N1(n_stocks,n_regions,n_ages);
  sim_log_N1.setZero();
  for(int s = 0; s < n_stocks; s++) if(N1_model(s) ==2)  for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,0)){
    Type mu = N1_repars(s,r,0);
    Type sigma = exp(N1_repars(s,r,1)); //marginal variance
    Type rho = geninvlogit(N1_repars(s,r,2),Type(-1),Type(1),Type(1));
    vector<Type> re_sr(n_ages);
    AR1(rho).simulate(re_sr);
    re_sr *= sigma;
    for(int a = 0; a < n_ages; a++) sim_log_N1(s,r,a) = re_sr(a) + mu;
  }
  return(sim_log_N1);
}

template <class Type>
vector<Type> get_SSB_y(int y, matrix<Type>NAA_spawn_y, array<Type> waa_ssb, array<Type> mature){
  /*
    provide annual SSB for each stock.
                    y: year index
          NAA_spawn_y: n_ages x n_stocks; numbers at age at time of spawning 
                  waa_ssb: (n_?) x n_years x n_ages_model. weight at age
               mature: n_stocks x n_years x n_ages; proportion mature
    
  */
  int n_stocks = NAA_spawn_y.cols();
  int n_ages = NAA_spawn_y.rows();
  
  vector<Type> SSB_y(n_stocks);// = get_SSB(NAA_ssb,waa_ssb,mature);
  SSB_y.setZero();
  // see(SSB_y);
  // see(n_stocks);
  // see(n_ages);
  // see(waa_ssb.dim);
  // see(mature.dim);
  // see(NAA_spawn_y.cols());
  // see(NAA_spawn_y.rows());
  // see(y);
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
    // see(NAA_spawn_y(a,s));
    // see(waa_ssb(s,y,a));
    // see(mature(s,y,a));
    // see(a);
    SSB_y(s) += NAA_spawn_y(a,s) * waa_ssb(s,y,a) * mature(s,y,a);
    // see(SSB_y);
  }
  // see(SSB_y);
  return SSB_y;
}

template <class Type>
matrix<Type> get_SSB(array<Type>NAA_spawn, array<Type> waa_ssb, array<Type> mature){
  /*
    provide annual SSB for each stock.
              NAA_spawn: n_stocks x n_years_pop x n_ages; numbers at age at time of spawning 
                  waa_ssb: (n_?) x n_years x n_ages_model. weight at age
               mature: n_stocks x n_years x n_ages; proportion mature
    
  */
  int n_stocks = NAA_spawn.dim(0);
  int n_y = NAA_spawn.dim(1);
  int n_ages = NAA_spawn.dim(2);
  
  matrix<Type> SSB(n_y, n_stocks);// = get_SSB(NAA_ssb,waa_ssb,mature);
  SSB.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_y; y++)
  {
    for(int a = 0; a < n_ages; a++) SSB(y,s) += NAA_spawn(s,y,a) * waa_ssb(s,y,a) * mature(s,y,a);
  }
  return(SSB);
}

template <class Type>
array<Type> get_NAA_1(vector<int> N1_model, array<Type> log_N1, array<int> NAA_where, array<Type> log_M, array<Type> FAA, 
  vector<int> which_F_age, vector<int> spawn_regions, 
  vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, array<Type> mu, 
  matrix<Type> L, vector<Type> fracyr_seasons, 
  int small_dim) {
  /* 
    get population age structure for the first year
             N1_model: (n_stocks) 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
               log_N1: (n_stocks x n_regions x n_ages) holding fixed or random effects paramters
            NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
           log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                  FAA: fishing mortality: n_fleets x n_years x n_ages
          which_F_age: (n_years_model + n_years_proj); which age of F to use for max F for Fmsy/Fxspr calculations and projections
        spawn_regions: n_stocks; which region spawning occurs for each stock
        fleet_regions: n_fleets; which region each fleet is operating
        fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
             can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                   mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                    L: n_years_model x n_regions; "extra" mortality rate
       fracyr_seasons: n_seasons: length of intervals for each season
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */
  int n_stocks = log_N1.dim(0);
  int n_fleets = FAA.dim(0);
  //int n_seasons = fleet_seasons.cols();
  int n_regions = log_N1.dim(1);
  int n_ages = log_M.dim(3);
  array<Type> NAA_1(n_stocks, n_regions, n_ages);
  NAA_1.setZero();
  //see("inside get_NAA_1");
  matrix<Type> sel1(n_fleets,n_ages);
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) sel1(f,a) = FAA(f,0,a);
  vector<Type> temp = sel1.colwise().sum();
  sel1 = sel1/temp(which_F_age(0)-1);
  //NAA_where(s,r,0) must be consistent with spawn_regions  
  for(int s = 0; s < n_stocks; s++) {
    if((N1_model(s) == 0) | (N1_model(s) == 2)) { //log_N1 is either fixed or random effects parameters for initial numbers at age
  //see("inside get_NAA_1 2");
      for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
        NAA_1(s,r,a) = exp(log_N1(s,r,a)); //log_N1 has to be mapped to not be estimated for NAA_where(s,r,a)==0
      }
    } else{ //N1_model(s) == 1
      array<Type> FAA1(n_fleets,1,n_ages);
      FAA1.setZero();
      for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
        FAA1(f,0,a) = exp(log_N1(s,spawn_regions(s)-1,1)) * sel1(f,a); //only 1 F0 per stock
      }
      array<Type> SAA1 = get_eq_SAA(0, fleet_regions, fleet_seasons, can_move, mig_type, FAA1, log_M, 
        mu, L, fracyr_seasons, small_dim);
      for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) if(NAA_where(s,i,a)) {
        NAA_1(s,i,a) += exp(log_N1(s,spawn_regions(s)-1,0)) * SAA1(s,a,spawn_regions(s)-1,i); //only 1 Rec per stock, this must be consistent with NAA_where
      }
    }
  }
  return NAA_1;
}

///////////////////THIS IS JUST TO return the components that go into creating equilibrium NAA for that option for initial numbers at age!!!!!!!!!!!!!!!
template <class Type>
vector< array<Type>> get_eq_NAA_components(vector<int> N1_model, array<Type> log_N1, array<int> NAA_where, array<Type> log_M, array<Type> FAA, 
  vector<int> which_F_age, vector<int> spawn_regions, 
  vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, array<Type> mu, 
  matrix<Type> L, vector<Type> fracyr_seasons, 
  int small_dim) {

  int n_stocks = log_N1.dim(0);
  int n_fleets = FAA.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_M.dim(3);
  array<Type> NAA_1(n_stocks, n_regions, n_ages);
  NAA_1.setZero();
  array<Type> sel1(n_fleets,n_ages);
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) sel1(f,a) = FAA(f,0,a);
  vector<Type> temp = sel1.matrix().colwise().sum();
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) sel1(f,a) = sel1(f,a)/temp(which_F_age(0)-1);
  vector< array<Type>> out(n_stocks*2+1);
  //NAA_where(s,r,0) must be consistent with spawn_regions  
  for(int s = 0; s < n_stocks; s++) {
    array<Type> FAA1(n_fleets,1,n_ages);
    FAA1.setZero();
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
      FAA1(f,0,a) = exp(log_N1(s,spawn_regions(s)-1,1)) * sel1(f,a); //only 1 F0 per stock
    }
    out(s*2) = FAA1;
    array<Type> SAA1 = get_eq_SAA(0, fleet_regions, fleet_seasons, can_move, mig_type, FAA1, log_M, 
      mu, L, fracyr_seasons, small_dim);
    out(s*2+1) = SAA1;
    for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) if(NAA_where(s,i,a)) {
      NAA_1(s,i,a) += exp(log_N1(s,spawn_regions(s)-1,0)) * SAA1(s,a,spawn_regions(s)-1,i); //only 1 Rec per stock, this must be consistent with NAA_where
    }
  }
  out(n_stocks*2) = sel1;
  return out;
}


template <class Type>
array<Type> get_NAA_y(int y, vector<int> NAA_re_model, array<Type> log_NAA, vector<int> N1_model, array<Type> log_N1, array<int> NAA_where, array<Type> log_M, 
  array<Type> FAA, vector<int> which_F_age, 
  vector<int> spawn_regions,
  vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, array<Type> mu, 
  array<Type> L, vector<Type> fracyr_seasons, int small_dim){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */
  /* 
    fill out numbers at age for year y
                    y: year index
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
              log_NAA: (n_stocks x n_regions x nyears-1 x n_ages) parameters for ages after year 1
             N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
               log_N1: (n_stocks x n_regions x n_ages) holding fixed or random effects paramters
            NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
           log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                  FAA: fishing mortality: n_fleets x n_years x n_seasons x n_ages
          which_F_age: (n_years_model + n_years_proj); which age of F to use for max F for Fmsy/Fxspr calculations and projections
        spawn_regions: n_stocks; which region spawning occurs for each stock
        fleet_regions: n_fleets; which region each fleet is operating
        fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
             can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                   mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                    L: n_years_model x n_regions; "extra" mortality rate
       fracyr_seasons: n_seasons: length of intervals for each season
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_y = log_M.dim(2); 
  int n_ages = log_M.dim(3);
  array<Type> NAA_y(n_stocks, n_regions, n_ages);
  NAA_y.setZero();
  if(y==0) {
    NAA_y = get_NAA_1(N1_model,log_N1, NAA_where, log_M, FAA, which_F_age, spawn_regions, fleet_regions, fleet_seasons, 
      can_move, mig_type, mu, L, fracyr_seasons, small_dim);
  } else{ 
    for(int s = 0; s < n_stocks; s++) {
      if(NAA_re_model(s) == 2){ //rec+1
        for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
          NAA_y(s,r,a) = exp(log_NAA(s,r,y-1,a));
        }
      }
      if(NAA_re_model(s) == 1) { //rec, rest of NAA will be 0. Need to populate with pred_NAA.
        NAA_y(s,spawn_regions(s)-1,0) = exp(log_NAA(s,spawn_regions(s)-1,y-1,0));
      }
    }
  }
  return(NAA_y);
}

template <class Type>
matrix<Type> get_NAA_spawn_y(int y, array<Type> NAA_y, array<Type> annual_SAA_spawn, vector<int> spawn_regions){
  int n_stocks = NAA_y.dim(0);
  int n_ages = NAA_y.dim(2);
  int n_regions = NAA_y.dim(1);
  matrix<Type> NAA_spawn_y(n_ages,n_stocks);
  NAA_spawn_y.setZero();
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++) {
    NAA_spawn_y(a,s) += NAA_y(s,r,a) * annual_SAA_spawn(s,y,a,r,spawn_regions(s)-1);
  }
  return NAA_spawn_y;
}

template <class Type>
array<Type> get_NAA_spawn(array<Type> NAA, array<Type> annual_SAA_spawn, vector<int> spawn_regions){
  int n_stocks = NAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);
  array<Type> NAA_spawn(n_stocks,n_years,n_ages);
  NAA_spawn.setZero();
  
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    NAA_spawn(s,y,a) += NAA(s,r,y,a) * annual_SAA_spawn(s,y,a,r,spawn_regions(s)-1);
  }
  return NAA_spawn;
}

template <class Type>
vector<Type> get_pred_recruit_y(int y, vector<int> recruit_model, matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, vector<int> NAA_re_model){
  /*
    provide "expected" recruitment (N(age 1)) for a given year
                  y: year (between 1 and n_years_model+n_years_proj)
      recruit_model: which recruitment model (1-4)
      mean_rec_pars: n_stocks x 2; recruitment parameters (defined in main code)
                SSB: n_years x n_stocks; of yearly SSB (uses y-1 for any S-R relationship)
                NAA: n_stocks x n_regions x n_years x n_ages; annual numbers at age by stock, region
           log_SR_a: yearly "a" parameters for SR function
           log_SR_b: yearly "b" parameters for SR function
         Ecov_how_R: integer vector with an element that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: region where spawning and recruitment occur for each stock
  */
  int n_stocks = NAA.dim(0);
  vector<Type> pred_recruit(n_stocks);
  pred_recruit.setZero();
  for(int s = 0; s < n_stocks; s++){
    if(NAA_re_model(s)>0) {//not SCAA
      if(recruit_model(s) == 1) { // random walk
        pred_recruit(s) = NAA(s,spawn_regions(s)-1,y-1,0);
      } else {
        if(recruit_model(s) == 2) {// random about mean
          pred_recruit(s) = exp(mean_rec_pars(s,0));
          int nE = Ecov_how_R.rows(); 
          for(int i=0; i < nE; i++){
            if(Ecov_how_R(i,s) == 1) pred_recruit(s) *= exp(Ecov_lm_R(s,y,i));
          }
        } else {
          if(recruit_model(s) == 3) {// BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
            pred_recruit(s) = exp(log_SR_a(y,s)) * SSB(y-1,s)/(1 + exp(log_SR_b(y,s))*SSB(y-1,s));
          } else {// recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
            pred_recruit(s) = exp(log_SR_a(y,s)) * SSB(y-1,s) * exp(-exp(log_SR_b(y,s)) * SSB(y-1,s));
          }
        }
      }
    }
  }
  return pred_recruit;
}

template <class Type>
vector<Type> get_pred_recruit_y(int y, vector<int> recruit_model, matrix<Type> mean_rec_pars, vector<Type> SSB_y_minus_1, 
  array<Type> NAA_y_minus_1, matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, vector<int> NAA_re_model){
  /*
    provide "expected" recruitment (N(age 1)) for a given year
                  y: year (between 1 and n_years_model+n_years_proj)
      recruit_model: which recruitment model (1-4)
      mean_rec_pars: recruitment parameters (defined in main code)
                SSB_y_minus: n_stocks; SSB (uses y-1 for any S-R relationship) at previous year
                NAA_y_minus_1: n_stocks x n_regions x n_ages; numbers at age by stock, region at previous year
           log_SR_a: yearly "a" parameters for SR function
           log_SR_b: yearly "b" parameters for SR function
         Ecov_how_R: integer vector with an element that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: region where spawning and recruitment occur for each stock
  */
  int n_stocks = NAA_y_minus_1.dim(0);
  vector<Type> pred_recruit(n_stocks);
  pred_recruit.setZero();
  for(int s = 0; s < n_stocks; s++){
    if(NAA_re_model(s)>0) {//not SCAA
      if(recruit_model(s) == 1) { // random walk
        pred_recruit(s) = NAA_y_minus_1(s,spawn_regions(s)-1,0);
      } else {
        if(recruit_model(s) == 2) {// random about mean
          pred_recruit(s) = exp(mean_rec_pars(s,0));
          int nE = Ecov_how_R.rows(); 
          for(int i=0; i < nE; i++){
            if(Ecov_how_R(i,s) == 1) pred_recruit(s) *= exp(Ecov_lm_R(s,y,i));
          }
        } else {
          if(recruit_model(s) == 3) {// BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
            pred_recruit(s) = exp(log_SR_a(y,s)) * SSB_y_minus_1(s)/(1 + exp(log_SR_b(y,s))*SSB_y_minus_1(s));
          } else {// recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
            pred_recruit(s) = exp(log_SR_a(y,s)) * SSB_y_minus_1(s) * exp(-exp(log_SR_b(y,s)) * SSB_y_minus_1(s));
          }
        }
      }
    }
  }
  return pred_recruit;
}

template <class Type>
array<Type> get_pred_N1(vector<int> N1_model, array<Type> N1, array<int> NAA_where, array<Type> N1_repars){
  /*
    provide the "expected" numbers at age in the first year. different from N1 only if N1 are random effects.
     N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
           N1: (n_stocks x n_regions x n_ages) parameters for estimated initial numbers at age
    NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
  */
  int n_stocks = N1.dim(0);
  int n_regions = N1.dim(1);
  int n_ages = N1.dim(2);
  array<Type> pred_N1(n_stocks, n_regions, n_ages);
  pred_N1.setZero();
  // see(n_stocks);
  // see(n_regions);
  // see(n_ages);
  // see(NAA_where.dim);
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
    if(N1_model(s) < 2) {
      pred_N1(s,r,a) = N1(s,r,a);
    } else {
      pred_N1(s,r,a) = exp(N1_repars(s,r,0)); //exp of mean of AR1 process ; bias correct?
    } 
  }
  return(pred_N1);
}

template <class Type>
array<Type> get_pred_NAA_y(int y, vector<int> N1_model, array<Type> N1, array<Type> N1_repars, array<int> NAA_where, vector<int> recruit_model, 
  matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> Ps, vector<int> NAA_re_model){

  /*
    provide "expected" numbers at age given NAA from previous time step (RECRUITMENT: ONLY FOR STOCKS with RE on NAA)
           N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
                 N1: n_stocks x n_regions x n_ages; array of parameters representing initial numbers at age.
          NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
      recruit_model: nstocks; which recruitment model (1-4)
      mean_rec_pars: n_stocks x 2; of any recruitment parameters (defined in main code)
                SSB: n_years x n_stocks; of yearly SSB (uses y-1 for any S-R relationship)
                NAA: nstocks x nregions x nyears x nages; array of numbers at age 
           log_SR_a: yearly "a" parameters for SR function for each stock
           log_SR_b: yearly "b" parameters for SR function for each stock
         Ecov_how_R: integer that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: nstocks; region where spawning and recruitment occur for each stock
                 Ps: array of annual PTMs
   */
  int n_stocks = NAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_ages = NAA.dim(3);
  array<Type> pred_NAA_y(n_stocks,n_regions,n_ages);
  pred_NAA_y.setZero();

  if(y==0) { //initial NAA
  //see("y=0");
    pred_NAA_y = get_pred_N1(N1_model, N1, NAA_where, N1_repars);
  } else {
  // Expected recruitment
    if(NAA_re_model.sum()>0){ //RE on NAA
      vector<Type> pred_recruit = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
        log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, NAA_re_model);
      for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0) pred_NAA_y(s,spawn_regions(s)-1,0) = pred_recruit(s);
    }
    // calculate pred_NAA for ages after recruitment
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
      for(int a = 1; a < n_ages; a++) {
        //pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
        pred_NAA_y(s,r,a) += Ps(s,y-1,a-1,rr,r) * NAA(s,rr,y-1,a-1);
      }
      //plus group
      //pred_NAA(n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
      pred_NAA_y(s,r,n_ages-1) += Ps(s,y-1,n_ages-1,rr,r) * NAA(s,rr,y-1,n_ages-1);
    }
  }
  return(pred_NAA_y);
}

template <class Type>
array<Type> get_pred_NAA_y(int y, vector<int> N1_model, array<Type> N1, array<Type> N1_repars, array<int> NAA_where, vector<int> recruit_model, 
  matrix<Type> mean_rec_pars, vector<Type> SSB_y_minus_1, array<Type> NAA_y_minus_1, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> Ps, vector<int> NAA_re_model){

  /*
    provide "expected" numbers at age given NAA from previous time step (RECRUITMENT: ONLY FOR STOCKS with RE on NAA)
           N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
                 N1: n_stocks x n_regions x n_ages; array of parameters representing initial numbers at age.
          NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
      recruit_model: nstocks; which recruitment model (1-4)
      mean_rec_pars: n_stocks x 2; of any recruitment parameters (defined in main code)
                SSB_y_minus_1: n_stocks; of yearly SSB at previous year
                NAA_y_minus_1: array of numbers at age at previous year
           log_SR_a: yearly "a" parameters for SR function for each stock
           log_SR_b: yearly "b" parameters for SR function for each stock
         Ecov_how_R: integer that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: nstocks; region where spawning and recruitment occur for each stock
                 Ps: array of annual PTMs
   */
  int n_stocks = NAA_y_minus_1.dim(0);
  int n_regions = NAA_y_minus_1.dim(1);
  int n_ages = NAA_y_minus_1.dim(2);
  array<Type> pred_NAA_y(n_stocks,n_regions,n_ages);
  pred_NAA_y.setZero();

  if(y==0) { //initial NAA
    pred_NAA_y = get_pred_N1(N1_model, N1, NAA_where, N1_repars);
  } else {
  // Expected recruitment ONLY FOR STOCKS with RE on NAA
    if(NAA_re_model.sum()>0){ //RE on NAA
      vector<Type> pred_recruit = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB_y_minus_1, NAA_y_minus_1, log_SR_a, 
      log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, NAA_re_model);
      for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0) pred_NAA_y(s,spawn_regions(s)-1,0) = pred_recruit(s);
    }
    // calculate pred_NAA for ages after recruitment
    for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++){
      for(int a = 1; a < n_ages; a++) {
        //pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
        pred_NAA_y(s,r,a) += Ps(s,y-1,a-1,rr,r) * NAA_y_minus_1(s,rr,a-1);
      }
      //plus group
      //pred_NAA(n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
      pred_NAA_y(s,r,n_ages-1) += Ps(s,y-1,n_ages-1,rr,r) * NAA_y_minus_1(s,rr,n_ages-1);
    }
  }
  return(pred_NAA_y);
}

template <class Type>
array<Type> get_pred_NAA(int N1_model, array<Type> N1, array<Type> N1_repars, array<int> NAA_where, vector<int> recruit_model, 
  matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> annual_Ps, int n_years_model, vector<int> NAA_re_model, matrix<Type> logR_proj){

  /*
    provide "expected" numbers at age given NAA from previous time step
           N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
                 N1: n_stocks x n_regions x n_ages; array of parameters representing initial numbers at age.
          NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
      recruit_model: nstocks; which recruitment model (1-4)
      mean_rec_pars: n_stocks x 2; of any recruitment parameters (defined in main code)
                SSB: n_years x n_stocks; of yearly SSB (uses y-1 for any S-R relationship)
                NAA: nstocks x nregions x nyears x nages; array of numbers at age 
           log_SR_a: yearly "a" parameters for SR function for each stock
           log_SR_b: yearly "b" parameters for SR function for each stock
         Ecov_how_R: integer that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: nstocks; region where spawning and recruitment occur for each stock
                 Ps: array of annual PTMs
   */
  int n_stocks = NAA.dim(0);
  int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);
  int n_regions = NAA.dim(1);
  array<Type> pred_NAA(n_stocks,n_regions,n_years,n_ages);
  pred_NAA.setZero();
  for(int y = 0; y < n_years; y++){
    array<Type> pred_NAA_y = get_pred_NAA_y(y, N1_model, N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB, NAA, 
      log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, NAA_re_model);
    for(int a = 0; a < n_ages; a++) for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
      //pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
      if((a==0) & (NAA_re_model(s)==0)) { //SCAA recruitment is not populated in get_pred_NAA_y
        if(r == spawn_regions(s)-1) {
          if(y > n_years_model - 1){ //projection year
            pred_NAA(s,r,y,a) = exp(logR_proj(y - n_years_model,s));
          }
          else pred_NAA(s,r,y,a) = NAA(s,r,y,a);
        }
      } else {
        pred_NAA(s,r,y,a) = pred_NAA_y(s,r,a);
      }
    }
  }
  return(pred_NAA);
}

template <class Type>
array<Type> get_all_NAA(vector<int> NAA_re_model, vector<int> N1_model, array<Type> N1, array<Type> N1_repars, 
  array<Type> log_NAA, array<int> NAA_where, 
  array<Type> mature, array<Type> waa_ssb,
  vector<int> recruit_model, matrix<Type> mean_rec_pars, matrix<Type> log_SR_a, matrix<Type> log_SR_b, 
  matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> annual_Ps, array<Type> annual_SAA_spawn, int n_years_model, int trace){
  /* 
    fill out numbers at age and "expected" numbers at age
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
             N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
               N1: (n_stocks x n_regions x n_ages) numbers at age in the first year
               N1:
               N1_repars:
              log_NAA: (n_stocks x n_regions x n_years_pop-1 x n_ages) parameters for ages after year 1
            NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
               mature: n_stocks x n_years_pop x n_ages; proportion mature
                  waa_ssb: (n_?) x n_years_pop x n_ages_model. weight at age
      recruit_model:
      mean_rec_pars:
      log_SR_a:
      log_SR_b:
      Ecov_how_R:
      Ecov_lm_R:
      spawn_regions:
      annual_Ps:
      annual_SAA_spawn:
      n_years_model: 
  */
  int n_stocks = log_NAA.dim(0);
  int n_regions = log_NAA.dim(1);
  int n_y = log_NAA.dim(2)+1; 
  int n_ages = log_NAA.dim(3);
  array<Type> NAA(2,n_stocks, n_regions, n_y, n_ages); //NAA AND pred_NAA
  NAA.setZero();
  array<Type> NAA_last = N1;
  if(trace) see("NAA1");
  matrix<Type> NAA_spawn_last = get_NAA_spawn_y(0, NAA_last, annual_SAA_spawn, spawn_regions);
  if(trace) see("NAA2");
  if(trace) see(NAA_spawn_last);
  vector<Type> SSB_last = get_SSB_y(0, NAA_spawn_last, waa_ssb, mature);
  if(trace) see("NAA3");
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
    NAA(0,s,r,0,a) = N1(s,r,a); //year 1 realized
  }
  if(trace){ 
    see("NAA4");
    see(N1_model);
    see(N1);
    see(N1_repars);
    see(NAA_where);
    see(recruit_model);
    see(mean_rec_pars);
    see(SSB_last);
    see(NAA_last);
    see(log_SR_a);
    see(log_SR_b);
    see(Ecov_how_R);
    see(Ecov_lm_R.dim);
    see(spawn_regions);
    see(annual_Ps.dim);
  }
  // vector<Type> logR_proj_y(n_stocks); //not used because no projection years
  // int is_projyr = 0;

  array<Type> pred_NAA_y = get_pred_NAA_y(0, N1_model, N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_last, NAA_last, 
    log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, NAA_re_model);
  if(trace) see("NAA5");
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
    NAA(1,s,r,0,a) = pred_NAA_y(s,r,a); //year 1 expected
  }
  if(trace) see("NAA6");

  for(int y = 1; y < n_years_model; y++){
    if(trace) see(y);
    pred_NAA_y = get_pred_NAA_y(y, N1_model, N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_last, NAA_last, 
      log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, NAA_re_model);
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
      if((a==0) & (NAA_re_model(s)==0)) { //SCAA recruitment is not populated in get_pred_NAA_y
        if(r == spawn_regions(s)-1) NAA(1,s,r,y,a) = exp(log_NAA(s,r,y-1,a));
      } else {
        NAA(1,s,r,y,a) = pred_NAA_y(s,r,a);
      }
      //NAA(1,s,r,y,a) = pred_NAA_y(s,r,a); //year y expected
    }
    if(trace) see("0.1");

    for(int s = 0; s < n_stocks; s++) {
      if(NAA_re_model(s) == 2){ //rec+1
        for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
          NAA(0,s,r,y,a) = exp(log_NAA(s,r,y-1,a)); //year y realized. rec+1
        }
      }
    if(trace) see("0.2");
      if(NAA_re_model(s) < 2) { //rec, Need to populate other ages with pred_NAA.
        //age 1 year y realized. rec
        NAA(0,s,spawn_regions(s)-1,y,0) = exp(log_NAA(s,spawn_regions(s)-1,y-1,0));
        //age 1 year y realized. SCAA
        //age 2+ year y realized. SCAA or rec
        for(int a = 1; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
          NAA(0,s,r,y,a) = pred_NAA_y(s,r,a);
        }
      }
    if(trace) see("0.3");
      for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) NAA_last(s,r,a) = NAA(0,s,r,y,a);
    if(trace) see("0.4");
      NAA_spawn_last = get_NAA_spawn_y(y, NAA_last,  annual_SAA_spawn, spawn_regions);
    if(trace) see("0.5");
    if(trace) see(NAA_spawn_last);
    if(trace) see(waa_ssb.dim);
    if(trace) see(mature.dim);
    if(trace) see(y);
    if(trace) see(SSB_last);
      SSB_last = get_SSB_y(y, NAA_spawn_last, waa_ssb, mature);
    if(trace) see("0.6");
    }
  }
  return(NAA);
}

template <class Type>
array<Type> update_all_NAA(int y, array<Type> all_NAA, vector<int> NAA_re_model, vector<int> N1_model, array<Type> N1, array<Type> N1_repars, 
  array<Type> log_NAA, array<int> NAA_where, 
  array<Type> mature, array<Type> waa_ssb,
  vector<int> recruit_model, matrix<Type> mean_rec_pars, matrix<Type> log_SR_a, matrix<Type> log_SR_b, 
  matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> annual_Ps, array<Type> annual_SAA_spawn, int n_years_model, matrix<Type> logR_proj, int proj_R_opt, matrix<Type> R_XSPR, 
  int bias_correct_pe, 
  array<Type> marg_NAA_sigma, 
  // array<Type> log_NAA_sigma, 
  int trace){
  /* 
    fill out numbers at age and "expected" numbers at age for year y (intended for projection years)
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
             N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
               N1: (n_stocks x n_regions x n_ages) numbers at age in the first year
               N1:
               N1_repars:
              log_NAA: (n_stocks x n_regions x n_years_pop-1 x n_ages) parameters for ages after year 1
            NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
               mature: n_stocks x n_years_pop x n_ages; proportion mature
                  waa: (n_?) x n_years_pop x n_ages_model. weight at age
      recruit_model:
      mean_rec_pars:
      log_SR_a:
      log_SR_b:
      Ecov_how_R:
      Ecov_lm_R:
      spawn_regions:
      annual_Ps:
      annual_SAA_spawn:
      n_years_model: 
  */
  if(trace) see(y);
  int n_stocks = log_NAA.dim(0);
  int n_regions = log_NAA.dim(1);
  int n_ages = log_NAA.dim(3);
  array<Type> updated_all_NAA = all_NAA;
  if(trace) see(updated_all_NAA.dim);
  array<Type> NAA_last(n_stocks,n_regions,n_ages);
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) NAA_last(s,r,a) = all_NAA(0,s,r,y-1,a);
  if(trace) see(NAA_last);
  
  matrix<Type> NAA_spawn_last = get_NAA_spawn_y(y-1, NAA_last,  annual_SAA_spawn, spawn_regions);
  if(trace) see(NAA_spawn_last);
  vector<Type> SSB_last = get_SSB_y(y-1, NAA_spawn_last, waa_ssb, mature);
  if(trace) see(SSB_last);

  array<Type> pred_NAA_y = get_pred_NAA_y(y, N1_model, N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_last, NAA_last, 
    log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, NAA_re_model);
  if(trace) see(pred_NAA_y);
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
    if((a==0) & (NAA_re_model(s)==0)) { //SCAA recruitment is not populated in get_pred_NAA_y
      if(r == spawn_regions(s)-1) pred_NAA_y(s,r,a) = exp(logR_proj(y-n_years_model,s)); // this function is called always in projection years
    } else {
      if((y>= n_years_model) & (proj_R_opt == 2)){ 
        //expected recruitment in projection years = RXSPR so that long term projections at FXSPR and SPR-based RFPs are consistent
        if((a == 0) & (r == spawn_regions(s)-1)) pred_NAA_y(s,r,a) = R_XSPR(y,s);
        if(bias_correct_pe) pred_NAA_y(s,r,a) *= exp(0.5 * pow(marg_NAA_sigma(s,r,a),2)); //take out bias correction in projections in this option
      }
    }
    updated_all_NAA(1,s,r,y,a) = pred_NAA_y(s,r,a);
  }
  if(trace) see("update_all_NAA(1)");

  for(int s = 0; s < n_stocks; s++) {
    if(NAA_re_model(s) == 2){ //rec+1
      for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
        updated_all_NAA(0,s,r,y,a) = exp(log_NAA(s,r,y-1,a)); //year y realized. rec+1
      }
    if(trace) see("NAA_re_model == 2, update_all_NAA(0)");
    }
    if(NAA_re_model(s) < 2) { //rec, Need to populate other ages with pred_NAA.
      //age 1 year y realized. rec
      if(NAA_re_model(s) == 1) { // projected recruitment is continued RE
        updated_all_NAA(0,s,spawn_regions(s)-1,y,0) = exp(log_NAA(s,spawn_regions(s)-1,y-1,0));
      } else { //SCAA
        //age 1 year y realized. SCAA
        updated_all_NAA(0,s,spawn_regions(s)-1,y,0) = exp(logR_proj(y-n_years_model,s));
      }
      //for SCAA or rec, age 2+ year y realized is deterministic
      for(int a = 1; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
        updated_all_NAA(0,s,r,y,a) = pred_NAA_y(s,r,a);
      }
    if(trace) see("NAA_re_model < 2, update_all_NAA(0)");
    }
  }
  return updated_all_NAA;
}

template <class Type>
array<Type> extract_NAA(array<Type> all_NAA){
  int n_stocks = all_NAA.dim(1);
  int n_regions = all_NAA.dim(2);
  int n_y = all_NAA.dim(3); 
  int n_ages = all_NAA.dim(4);
  array<Type> NAA(n_stocks,n_regions,n_y,n_ages);
  NAA.setZero();
  for(int y = 0; y < n_y; y++) for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    NAA(s,r,y,a) = all_NAA(0,s,r,y,a);
  }
  return NAA;
}

template <class Type>
array<Type> extract_pred_NAA(array<Type> all_NAA){
  int n_stocks = all_NAA.dim(1);
  int n_regions = all_NAA.dim(2);
  int n_y = all_NAA.dim(3); 
  int n_ages = all_NAA.dim(4);
  array<Type> pred_NAA(n_stocks,n_regions,n_y,n_ages);
  pred_NAA.setZero();
  for(int y = 0; y < n_y; y++) for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    pred_NAA(s,r,y,a) = all_NAA(1,s,r,y,a);
  }
  return pred_NAA;
}


template <class Type>
array<Type> get_NAA_devs(array<Type> all_NAA, array<int> NAA_where, vector<int> NAA_re_model){
  int n_stocks = all_NAA.dim(1);
  int n_regions = all_NAA.dim(2);
  int n_y = all_NAA.dim(3); 
  int n_ages = all_NAA.dim(4);
  array<Type> NAA_devs(n_stocks,n_regions,n_y,n_ages);
  NAA_devs.setZero();
  // see("get_NAA_devs");
  // see(all_NAA.dim);
  // see(n_stocks);
  // see(n_regions);
  // see(n_y);
  // see(n_ages);
  for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0) {
    for(int r = 0; r < n_regions; r++) for(int y = 1; y < n_y; y++) for(int a = 0; a < n_ages; a++) {
      // see(y);
      // see(s);
      // see(a);
      // see(r);
      // see(all_NAA(0,s,r,y,a));
      // see(all_NAA(1,s,r,y,a));
      if(NAA_where(s,r,a)) NAA_devs(s,r,y,a) = log(all_NAA(0,s,r,y,a)) - log(all_NAA(1,s,r,y,a));
    }
  }
  // see("end get_NAA_devs");
  return NAA_devs;
}

template <class Type>
array<Type> get_NAA_y(int y, array<Type> NAA){
  int n_stocks = NAA.dim(0);
  int n_regions = NAA.dim(1);
  int n_ages = NAA.dim(3);
  array<Type> NAA_y(n_stocks,n_regions,n_ages);
  NAA_y.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    NAA_y(s,r,a) = NAA(s,r,y,a);
  }
  return NAA_y;
}

template <class Type>
array<Type> get_log_NAA_rep(array<Type> NAA, array<int> NAA_where){
  array<Type> log_NAA = NAA;
  log_NAA.setZero();
  for(int s = 0; s < NAA.dim(0); s++) for(int r = 0; r < NAA.dim(1); r++) for(int y = 0; y < NAA.dim(2); y++) for(int a = 0; a < NAA.dim(3); a++){
    if(NAA_where(s,r,a)) {
      log_NAA(s,r,y,a) = log(NAA(s,r,y,a));
    }
    else log_NAA(s,r,y,a) = -100.0; //no NAA there
  }
  return log_NAA;
}

template <class Type>
matrix<Type> get_SR_log_a(vector<int> recruit_model, matrix<Type> mean_rec_pars, array<Type> Ecov_lm_R, matrix<int> Ecov_how_R){
  /*
    make annual stock recruit log(a) parameters for each stock
      recruit_model: n_stocks; which recruitment model; 3=BH, 4=Ricker
      mean_rec_pars: n_stocks x 2; base recruitment parameters
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
         Ecov_how_R: n_Ecov x n_stocks: specific to recruitment effects. 0 = no effect, 1 = controlling, 2 = limiting, 3 = lethal, 4 = masking, 5 = directive
  */

  int n_y = Ecov_lm_R.dim(1);  
  int n_s = Ecov_lm_R.dim(0);
  int n_Ecov = Ecov_lm_R.dim(2);
  matrix<Type> log_SR_a(n_y,n_s);
  log_SR_a.setZero();
  for(int s = 0; s < n_s; s++) {
    for(int y = 0; y < n_y; y++) {
      log_SR_a(y,s) += mean_rec_pars(s,0);
      if(recruit_model(s) == 3) {//BH stock recruit
        for(int i = 0; i < n_Ecov; i++) {
          // (1) "controlling" = dens-indep mortality or (4) "masking" = metabolic/growth (decreases dR/dS)
          if((Ecov_how_R(i,s) == 1) | (Ecov_how_R(i,s) == 4)) {
            log_SR_a(y,s) += Ecov_lm_R(s,y,i); //only first "region" used
          }
        }
      }
      if(recruit_model(s) == 4) { //Ricker stock recruit
        for(int i = 0; i < n_Ecov; i++) {
          if(Ecov_how_R(i,s) == 1) { // "controlling" = dens-indep mortality         
            log_SR_a(y,s) += Ecov_lm_R(s,y,i);
          }
        }
      }
    }
  }
  return log_SR_a;
}

template <class Type>
matrix<Type> get_SR_log_b(vector<int> recruit_model, matrix<Type> mean_rec_pars, array<Type> Ecov_lm_R, matrix<int> Ecov_how_R){
  /*
    make annual stock recruit log(b) parameters for each stock
      recruit_model: n_stocks; which recruitment model; 3=BH, 4=Ricker
      mean_rec_pars: n_stocks x 2; base recruitment parameters
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
         Ecov_how_R: n_Ecov x n_stocks: specific to recruitment effects. 0 = no effect, 1 = controlling, 2 = limiting, 3 = lethal, 4 = masking, 5 = directive
  */
  int n_y = Ecov_lm_R.dim(1);  
  int n_s = Ecov_lm_R.dim(0);
  int n_Ecov = Ecov_lm_R.dim(2);
  matrix<Type> log_SR_b(n_y,n_s);
  log_SR_b.setZero();
  for(int s = 0; s < n_s; s++) {
    for(int y = 0; y < n_y; y++) {
      log_SR_b(y,s) += mean_rec_pars(s,1);
      if(recruit_model(s) == 3) { //BH stock recruit
        for(int i = 0; i < n_Ecov; i++) {
          // (1) "controlling" = dens-indep mortality or (4) "masking" = metabolic/growth (decreases dR/dS)
          // (2) "limiting" = carrying capacity or (4) "masking" = metabolic/growth (decreases dR/dS)
          if((Ecov_how_R(i,s) == 2) | (Ecov_how_R(i,s) == 4)) {
            log_SR_b(y,s) += Ecov_lm_R(s,y,i);
          }
        }
      }
      if(recruit_model(s) == 4) //Ricker stock recruit
      {
        for(int i = 0; i < n_Ecov; i++) {
          if(Ecov_how_R(i,s) == 4) { // "masking" = metabolic/growth (decreases dR/dS)
            //NB: this is not identical to Iles and Beverton (1998), but their definition can give negative values of "b"
            log_SR_b(y,s) += 1.0 + Ecov_lm_R(s,y,i);
          }
        }
      }
    }
  }
  return log_SR_b;
}

template <class Type>
array<Type> get_NAA_index(array<Type> NAA, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, 
  vector<Type> fracyr_seasons,
  matrix<Type> fracyr_indices, vector<int> index_seasons, vector<int> index_regions, array<Type> FAA, array<Type> log_M, 
  array<Type> mu, matrix<Type> L, int n_years_model){
  /*
    produce the annual survival probabilities up to time of spawning for a given stock, age, season, year
                NAA: nstocks x nregions x nyears x nages; array of numbers at age 
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_seasons: n_seasons; length of intervals for each season
      fracyr_indices: n_indices; length of intervals for each index
      index_seasons: n_indices; which season the index occurs in
      index_regions: n_indices: which region the index is observing
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_indices = index_seasons.size();
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  //int n_years = log_M.dim(2);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim

  array<Type> NAA_index(n_stocks,n_indices,n_years_model,n_ages);
  NAA_index.setZero();
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) {
      for(int i = 0; i < n_indices; i++) {
        if(t == index_seasons(i)-1){ 
          //P(0,t) x P(t_i-t): PTM over interval from beginning of year to time of index within the season
          matrix<Type> P_index = P_y * get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_indices(y,i), 
            FAA, log_M, mu, L);
          for(int r = 0; r < n_regions; r++) NAA_index(s,i,y,a) += P_index(r,index_regions(i)-1) * NAA(s,r,y,a);
        }
      }
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
  }
  return(NAA_index);
}

template <class Type>
array<Type> get_NAA_catch(array<Type> NAA, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, 
  vector<Type> fracyr_seasons, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the numbers caught by stock, fleet, year, season, age up to time of spawning for a given stock, age, season, year
                NAA: nstocks x nregions x nyears x nages; array of numbers at age 
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons: 0/1 indicator whether fleet is operating in a given season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_seasons: n_seasons; length of intervals for each season
      fracyr_indices: n_indices; length of intervals for each index
      index_seasons: n_indices; which season the index occurs in
      index_regions: n_indices: which region the index is observing
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_years = log_M.dim(2);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim

  array<Type> NAA_catch(n_stocks,n_fleets,n_years,n_seasons,n_ages);
  NAA_catch.setZero();
  //array<Type> annual_Ps_index(n_stocks,n_years,n_ages,P_dim,P_dim);
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();  
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) {
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      if(sum(vector<int> (fleet_seasons.col(t)))>0){
        array<Type> NAA_alive(n_stocks,n_regions,n_ages);
        NAA_alive.setZero();
        //number alive a the beginning of this season
        for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
          NAA_alive(s,rr,a) += P_y(r,rr) * NAA(s,r,y,a);
        }
        for(int f = 0; f < n_fleets; f++) if(fleet_seasons(f,t)){
          //number caught during this season
          for(int r = 0; r < n_regions; r++){
            NAA_catch(s,f,y,t,a) += NAA_alive(s,r,a) * P_t(r,n_regions + f);
          }
        }
      }
      P_y = P_y * P_t;
    }
  }
  return(NAA_catch);
}


template <class Type>
matrix<Type> get_NAA_nll(vector<int> NAA_re_model, array<Type> all_NAA, array<Type> log_NAA_sigma, array<Type> trans_NAA_rho, 
  array<int> NAA_where,
  vector<int> spawn_regions, vector<int> years_use, int bias_correct_pe, int decouple_recruitment = 0, int use_alt_AR1 = 0){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */
  //currently independent 2D (at most) AR1 processes by stock and region. not all ages may be available in all regions. number of stocks will typically be small
  //NAA_logsigma n_stocks x n_ages x n_regions
  //trans_NAA_rho n_stocks x n_regions x 3 (rho_a, rho_y, recruits rho_y) 
  //years_use is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  // array<Type> all_NAA = get_all_NAA(NAA_re_model, N1_model, N1, N1_repars, log_NAA, NAA_where, 
  //   mature, waa_ssb, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
  //   Ecov_how_R, Ecov_lm_R, spawn_regions,  annual_Ps, annual_SAA_spawn, n_years_model,0); //log_NAA should be mapped accordingly to exclude NAA=0 e.g., recruitment by region.
  array<Type> NAA = extract_NAA(all_NAA);
  array<Type> pred_NAA = extract_pred_NAA(all_NAA);

  int n_stocks = NAA.dim(0);
  int n_years = years_use.size();
  //int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);
  int n_regions = NAA.dim(1);
  int rho_y_ind = 1;
  if(decouple_recruitment) rho_y_ind = 2;

  vector<Type> marginal_sigma(n_ages); //sigmas for one stock
  Type NAA_rho_a = 0, NAA_rho_y = 0;
  matrix<Type> nll_NAA(n_stocks,n_regions);
  nll_NAA.setZero();
  //can we do vector< vector< matrix<Type>>>?
  //NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0){
    marginal_sigma.setZero(); //clear for each stock
    if((NAA_re_model(s) == 1) | ((NAA_re_model(s) == 2) & decouple_recruitment)){ //"rec"
      vector<Type> NAA_devs_r_s(n_years-1);
      // for NAA_re_model = 1, must make sure that rho_a = 0 and rho_y is set appropriately (cor = "iid" or "ar1_y") on R side
      NAA_rho_y = geninvlogit(trans_NAA_rho(s,spawn_regions(s)-1,rho_y_ind), Type(-1), Type(1), Type(1)); //using scale =1 ,2 is legacy
      marginal_sigma(0) = exp(log_NAA_sigma(s,spawn_regions(s)-1,0)) * pow(1-pow(NAA_rho_y,2),-0.5);
      for(int y = 1; y < n_years; y++) NAA_devs_r_s(y-1) = log(NAA(s,spawn_regions(s)-1,years_use(y),0)) - log(pred_NAA(s,spawn_regions(s)-1,years_use(y),0));
      if(bias_correct_pe) NAA_devs_r_s += 0.5*pow(marginal_sigma(0),2); //make sure this is ok when just recruitment is random.
      if(use_alt_AR1==0){
        nll_NAA(s,spawn_regions(s)-1) += SCALE(AR1(NAA_rho_y), marginal_sigma(0))(NAA_devs_r_s);
      } else {
        nll_NAA(s,spawn_regions(s)-1) +=  dar1(NAA_devs_r_s, trans_NAA_rho(s,spawn_regions(s)-1,rho_y_ind), log_NAA_sigma(s,spawn_regions(s)-1,0), 0);
      }
    }
    if(NAA_re_model(s) == 2){ //"rec+1"
      for(int r = 0; r < n_regions; r++){

        // for NAA_re_model = 1, must make sure that rho_a = 0 and rho_y is set appropriately (cor = "iid" or "ar1_y") on R side
        NAA_rho_a = geninvlogit(trans_NAA_rho(s,r,0), Type(-1), Type(1), Type(1)); //using scale =1 ,2 is legacy
        NAA_rho_y = geninvlogit(trans_NAA_rho(s,r,1), Type(-1), Type(1), Type(1)); //using scale =1 ,2 is legacy
        int n_age_s_r = 0;
        int age_start = 0;
        if(decouple_recruitment) age_start = 1;
        for(int a = age_start; a< n_ages; a++) {
          if(NAA_where(s,r,a)) n_age_s_r++; //start after recruitment when decouple_recruitment= 1
        }
        if(n_age_s_r>0) {// has to be some fish of some age in this region
          array<Type> NAA_devs_s_r(n_years-1, n_age_s_r);
          NAA_devs_s_r.setZero();
          vector<Type> marginal_sigma_s_r(n_age_s_r), log_sigma_s_r(n_age_s_r);
          int k=0;
          for(int a = age_start; a< n_ages; a++) if(NAA_where(s,r,a)) {
            log_sigma_s_r(k) = log_NAA_sigma(s,r,a);
            marginal_sigma_s_r(k) = exp(log_sigma_s_r(k)) * pow((1-pow(NAA_rho_y,2))*(1-pow(NAA_rho_a,2)),-0.5);
            k++;
          }
          for(int y = 1; y < n_years; y++) {
            k=0;
            for(int a = age_start; a< n_ages; a++) if(NAA_where(s,r,a)) {
              NAA_devs_s_r(y-1,k) = log(NAA(s,r,years_use(y),a)) - log(pred_NAA(s,r,years_use(y),a));
              if(bias_correct_pe) NAA_devs_s_r(y-1,k) += 0.5*pow(marginal_sigma_s_r(k),2);
              k++;
            }
          }
          // see(NAA_devs_s_r);
          // see(NAA_devs_s_r.rows());
          // see(NAA_devs_s_r.cols());
          // see(marginal_sigma_s_r);
          // see(NAA_rho_a);
          // see(NAA_rho_y);
          if(use_alt_AR1==0){
            nll_NAA(s,r) += SEPARABLE(VECSCALE(AR1(NAA_rho_a), marginal_sigma_s_r),AR1(NAA_rho_y))(NAA_devs_s_r);
          } else {
            // see("using alt2dar1");
            nll_NAA(s,r) += d2dar1(NAA_devs_s_r, trans_NAA_rho(s,r,1), trans_NAA_rho(s,r,0), log_sigma_s_r, 0);
          }
          // see(nll_NAA(s,r));
          //std::exit(EXIT_FAILURE);
        }
      }
    }
  }
  return nll_NAA;
}

template <class Type>
array<Type> simulate_NAA_devs(array<Type> NAA_devs, vector<int> NAA_re_model, array<Type> log_NAA_sigma, array<Type> trans_NAA_rho, array<int> NAA_where, 
  vector<int> spawn_regions, vector<int> years_use, int bias_correct_pe, int decouple_recruitment = 0, int use_alt_AR1 = 0, int ystart = 0){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */
  //simulate NAA devs for all model and projection years. The predicted NAA in projection years can change and simulate_log_NAA will change in projection years
  //currently independent 2D (at most) AR1 processes by stock and region. not all ages may be available in all regions. number of stocks will typically be small
  //NAA_logsigma n_stocks x n_ages x n_regions
  //trans_NAA_rho n_stocks x n_regions x 3 (rho_a, rho_y, recruit rho_y) 
  //years_use is possibly a subset of years to use for evaluating likelihood (and simulating values). normally = 0,....,n_years_model-1
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  int n_stocks = NAA_devs.dim(0);
  //int n_years_sim = years_use.size()-1;
  int n_years = years_use.size();
  //int n_years_pop = NAA_devs.dim(2); //all years
  int n_ages = NAA_devs.dim(3);
  int n_regions = NAA_devs.dim(1);
  Type NAA_rho_y = 0, NAA_rho_a = 0;
  array<Type> NAA_devs_out = NAA_devs; //(n_stocks, n_regions, n_years_pop, n_ages); //same dims as that provided by get_NAA_devs
  vector<Type> marginal_sigma(n_ages);

  int rho_y_ind = 1;
  if(decouple_recruitment) rho_y_ind = 2;

  //NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0){
    marginal_sigma.setZero();  //clear for each stock
    if((NAA_re_model(s) == 1) | ((NAA_re_model(s) == 2) & decouple_recruitment)){ //"rec"
      // see(decouple_recruitment);
      // see(NAA_re_model(s));
      // for NAA_re_model = 1, must make sure that rho_a = 0 and rho_y is set appropriately (cor = "iid" or "ar1_y") on R side
      NAA_rho_y = geninvlogit(trans_NAA_rho(s,spawn_regions(s)-1,rho_y_ind), Type(-1), Type(1), Type(1)); //using scale =1 ,2 is legacy
      marginal_sigma(0) = exp(log_NAA_sigma(s,spawn_regions(s)-1,0)) * pow(1-pow(NAA_rho_y,2),-0.5);
      vector<Type> NAA_devs_r_s(n_years-1);
      NAA_devs_r_s.setZero();
      if(use_alt_AR1==1){ //do simulation using conditional pdfs
        for(int y = 1; y < n_years; y++){
          NAA_devs_r_s(y-1) = NAA_devs_out(s,spawn_regions(s)-1,years_use(y),0);
        }
        NAA_devs_r_s = rar1(NAA_devs_r_s,trans_NAA_rho(s,spawn_regions(s)-1,rho_y_ind),log_NAA_sigma(s,spawn_regions(s)-1,0),0,ystart,bias_correct_pe);
      } else{
        AR1(NAA_rho_y).simulate(NAA_devs_r_s); // sigma = 1, scale below
        NAA_devs_r_s = marginal_sigma(0) * NAA_devs_r_s;
        if(bias_correct_pe) NAA_devs_r_s -= 0.5*pow(marginal_sigma(0),2);
      }
      for(int y = ystart+1; y < n_years; y++){
        NAA_devs_out(s,spawn_regions(s)-1,years_use(y),0) = NAA_devs_r_s(y-1);
      }
    }
    if(NAA_re_model(s) == 2){ //"rec+1"
      for(int r = 0; r < n_regions; r++){
        NAA_rho_a = geninvlogit(trans_NAA_rho(s,r,0), Type(-1), Type(1), Type(1)); //using scale =1 ,2 is legacy
        NAA_rho_y = geninvlogit(trans_NAA_rho(s,r,1), Type(-1), Type(1), Type(1)); //using scale =1 ,2 is legacy
        int n_age_s_r = 0;
        int age_start = 0;
        if(decouple_recruitment) age_start = 1;
        for(int a = age_start; a< n_ages; a++) {
          if(NAA_where(s,r,a)) n_age_s_r++; //start after recruitment when decouple_recruitment= 1
        }
        if(n_age_s_r>0) {// has to be some fish of some age in this region
          array<Type> NAA_devs_s_r(n_years-1, n_age_s_r);
          NAA_devs_s_r.setZero();
          vector<Type> marginal_sigma_s_r(n_age_s_r), log_sigma_s_r(n_age_s_r);
          int k=0;
          for(int a = age_start; a< n_ages; a++) if(NAA_where(s,r,a)) {
            log_sigma_s_r(k) = log_NAA_sigma(s,r,a);
            marginal_sigma_s_r(k) = exp(log_sigma_s_r(k)) * pow((1-pow(NAA_rho_y,2))*(1-pow(NAA_rho_a,2)),-0.5);
            k++;
          }
          // see(use_alt_AR1);
          // see(ystart);
          // see(age_start);
          if(use_alt_AR1==1){ //do simulation using conditional pdfs
            for(int y = ystart+1; y < n_years; y++) {
              k=0;
              for(int a = age_start; a< n_ages; a++) if(NAA_where(s,r,a)) {
                // if(a == 0) see(NAA_devs_out(s,r,years_use(y),a));
                NAA_devs_s_r(y-1,k) = NAA_devs_out(s,r,years_use(y),a);
                k++;
              }
              NAA_devs_s_r = r2dar1(NAA_devs_s_r,trans_NAA_rho(s,r,1), trans_NAA_rho(s,r,0), log_sigma_s_r,0,ystart,bias_correct_pe);
            }
          } else {
            SEPARABLE(VECSCALE(AR1(NAA_rho_a), marginal_sigma_s_r),AR1(NAA_rho_y)).simulate(NAA_devs_s_r); // scaled here
            for(int y = ystart+1; y < n_years; y++) { 
              k=0;
              for(int a = age_start; a< n_ages; a++) if(NAA_where(s,r,a)) {
                if(bias_correct_pe) NAA_devs_s_r(y-1,k) -= 0.5*pow(marginal_sigma_s_r(k),2);
                k++;
              }
            }
          }
          for(int y = ystart+1; y < n_years; y++) { 
            k=0;
            for(int a = age_start; a< n_ages; a++) if(NAA_where(s,r,a)) {
              // if(a == 0) see(years_use(y));
              // if(a == 0) see(NAA_devs_out(s,r,years_use(y),a));
              NAA_devs_out(s,r,years_use(y),a) = NAA_devs_s_r(y-1,k);
              // if(a == 0) see(NAA_devs_out(s,r,years_use(y),a));
              k++;
            }
          }
        }
      }
    }
  }
  return NAA_devs_out;
}

// stopped here
template <class Type>
matrix<Type> get_simulated_log_NAA(vector<int> N1_model, array<Type> N1, array<Type> N1_repars, vector<int> NAA_re_model, array<Type> NAA_devs, array<Type> log_NAA,
  array<int> NAA_where, vector<int> recruit_model, matrix<Type> mean_rec_pars, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R,
  vector<int> spawn_regions, array<Type> annual_Ps, array<Type> annual_SAA_spawn, array<Type> waa_ssb, 
  array<Type> mature, int n_years_model, matrix<Type> logR_proj){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */
  //does not simulate anything, but uses simulated NAA_devs to construct log_NAA, needs to be iteratively called in projection years given 
  // updates in FAA, MAA, annual_Ps, etc.
  // only needed for NAA_re_model = 1 or 2
  // NAA_devs is n_stocks x n_regions x n_years_pop x n_ages(ish). Values in first year are 0. 
  //currently independent 2D (at most) AR1 processes by stock and region. not all ages may be available in all regions. number of stocks will 
  // typically be small.
  //NAA_logsigma n_stocks x n_ages x n_regions
  //NAA_trans_rho n_stocks x 4 (rho_a, rho_y, rho_r, rho_s) 

  int n_stocks = log_NAA.dim(0);
  int n_regions = log_NAA.dim(1);
  int n_years_pop = log_NAA.dim(2)+1;
  int n_ages = log_NAA.dim(3);
  array<Type> sim_log_NAA = log_NAA;
  sim_log_NAA.setZero();
  array<Type> NAA_y_minus_1 = N1;
  matrix<Type> NAA_spawn_y_minus_1 = get_NAA_spawn_y(0, NAA_y_minus_1,  annual_SAA_spawn, spawn_regions);
  vector<Type> SSB_y_minus_1 = get_SSB_y(0, NAA_spawn_y_minus_1, waa_ssb, mature);
  for(int y = 1; y < n_years_pop; y++) {

    array<Type> pred_NAA_y = get_pred_NAA_y(y, N1_model, N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_y_minus_1, NAA_y_minus_1, 
      log_SR_a, log_SR_b, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps, NAA_re_model);
    NAA_y_minus_1.setZero();
    for(int s = 0; s < n_stocks; s++){
      if(NAA_re_model(s)==1){ //rec
        sim_log_NAA(s,spawn_regions(s)-1,y-1,0) += log(pred_NAA_y(s,spawn_regions(s)-1,0)) + NAA_devs(s,spawn_regions(s)-1,y,0);
        NAA_y_minus_1(s,spawn_regions(s)-1,0) = exp(sim_log_NAA(s,spawn_regions(s)-1,y-1,0));
        for(int a = 1; a < n_ages; a++) for(int r = 0; r < n_regions; r++) NAA_y_minus_1(s,r,a) = pred_NAA_y(s,r,a);
      }
      if(NAA_re_model(s)==2){ //rec+1
        for(int r = 0; r < n_regions; r++){
          for(int a = 0; a < n_ages; a++) if(NAA_where(s,r,a)) {
            sim_log_NAA(s,r,y-1,a) += log(pred_NAA_y(s,r,a)) + NAA_devs(s,r,y,a);
            NAA_y_minus_1(s,r,a) = exp(sim_log_NAA(s,r,y-1,a));
          }
        }
      }
    }
    NAA_spawn_y_minus_1 = get_NAA_spawn_y(y, NAA_y_minus_1, annual_SAA_spawn, spawn_regions);
    SSB_y_minus_1 = get_SSB_y(y, NAA_spawn_y_minus_1, waa_ssb, mature);
  }

  return sim_log_NAA;
}
