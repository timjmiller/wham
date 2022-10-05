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
vector<Type> get_SSB_y(int y, matrix<Type>NAA_spawn_y, array<Type> waa, vector<int> waa_pointer_ssb, array<Type> mature){
  /*
    provide annual SSB for each stock.
                    y: year index
          NAA_spawn_y: n_ages x n_stocks; numbers at age at time of spawning 
                  waa: (n_?) x n_years x n_ages_model. weight at age
      waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
               mature: n_stocks x n_years x n_ages; proportion mature
    
  */
  int n_stocks = NAA_spawn_y.cols();
  int n_ages = NAA_spawn_y.rows();
  
  vector<Type> SSB_y(n_stocks);// = get_SSB(NAA_ssb,waa,waa_pointer_ssb,mature);
  SSB_y.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
    SSB_y(s) += NAA_spawn_y(a,s) * waa(waa_pointer_ssb(s)-1,y,a) * mature(s,y,a);
  }
  return(SSB_y);
}

template <class Type>
matrix<Type> get_SSB(array<Type>NAA_spawn, array<Type> waa, vector<int> waa_pointer_ssb, array<Type> mature){
  /*
    provide annual SSB for each stock.
              NAA_spawn: n_stocks x n_years_pop x n_ages; numbers at age at time of spawning 
                  waa: (n_?) x n_years x n_ages_model. weight at age
      waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
               mature: n_stocks x n_years x n_ages; proportion mature
    
  */
  int n_stocks = NAA_spawn.dim(0);
  int n_y = NAA_spawn.dim(1);
  int n_ages = NAA_spawn.dim(2);
  
  matrix<Type> SSB(n_y, n_stocks);// = get_SSB(NAA_ssb,waa,waa_pointer_ssb,mature);
  SSB.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_y; y++)
  {
    for(int a = 0; a < n_ages; a++) SSB(y,s) += NAA_spawn(s,y,a) * waa(waa_pointer_ssb(s)-1,y,a) * mature(s,y,a);
  }
  return(SSB);
}

template <class Type>
array<Type> get_NAA_1(vector<int> N1_model, array<Type> log_N1, array<int> NAA_where, array<Type> log_M_base, array<Type> FAA, 
  matrix<int> which_F_age, vector<int> spawn_seasons, 
  vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, matrix<Type> fracyr_SSB,  array<Type> mu, 
  matrix<Type> L,  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, 
  vector<int> avg_years_ind, int small_dim) {
  /* 
    get population age structure for the first year
             N1_model: (n_stocks) 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
               log_N1: (n_stocks x n_regions x n_ages) holding fixed or random effects paramters
            NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
           log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                  FAA: fishing mortality: n_fleets x n_years x n_ages
          which_F_age: (n_years_model + n_years_proj x 2); which age,fleet of F to use for max F for msy/ypr calculations and projections
        spawn_seasons: n_stocks; which season spawning occurs for each stock
        fleet_regions: n_fleets; which region each fleet is operating
        fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
             can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
           fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                   mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                    L: n_years_model x n_regions; "extra" mortality rate
               mature: n_stocks x n_years x n_ages; proportion mature
                  waa: (n_?) x n_years x n_ages_model. weight at age
      waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
       fracyr_seasons: n_seasons: length of intervals for each season
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */
  int n_stocks = log_N1.dim(0);
  int n_fleets = FAA.dim(0);
  int n_seasons = fleet_seasons.cols();
  int n_regions = log_N1.dim(1);
  int n_ages = log_M_base.dim(3);
  array<Type> NAA_1(n_stocks, n_regions, n_ages);
  NAA_1.setZero();
  matrix<Type> sel1 = get_sel_proj(0, FAA, avg_years_ind, which_F_age);

  //NAA_where(s,r,0) must be consistent with spawn_regions  
  for(int s = 0; s < n_stocks; s++) {
    if((N1_model(s) == 0) | (N1_model(s) == 2)) { //log_N1 is either fixed or random effects parameters for initial numbers at age
      for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
        NAA_1(s,r,a) = exp(log_N1(s,r,a)); //log_N1 has to be mapped to not be estimated for NAA_where(s,r,a)==0
      }
    } else{ //N1_model(s) == 2
      vector<int> no_avg_yrs_ind(1);
      no_avg_yrs_ind(0) = 0;
      array<Type> FAA1(n_fleets,1,n_ages);
      FAA1.setZero();
      for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
        FAA1(f,0,a) = exp(log_N1(s,0,1)) * sel1(f,a); //only 1 F0 per stock
      }
      array<Type> SAA1 = get_eq_SAA(0, spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB, FAA1, log_M_base, 
        mu, L, mature, waa, waa_pointer_ssb, fracyr_seasons, small_dim);
      for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) if(NAA_where(s,i,a)) {
        NAA_1(s,i,a) += exp(log_N1(s,i,0)) * SAA1(s,i,i,a); //only 1 Rec per stock, this must be consistent with NAA_where
      }
    }
  }
  return NAA_1;
}

template <class Type>
array<Type> get_NAA_y(int y, vector<int> NAA_re_model, array<Type> log_NAA, vector<int> N1_model, array<Type> log_N1, array<int> NAA_where, array<Type> log_M_base, 
  array<Type> FAA, matrix<int> which_F_age, vector<int> spawn_seasons, 
  vector<int> spawn_regions,
  vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, matrix<Type> fracyr_SSB,  array<Type> mu, 
  array<Type> L, array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, 
  vector<int> avg_years_ind, int small_dim){
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
           log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                  FAA: fishing mortality: n_fleets x n_years x n_seasons x n_ages
          which_F_age: (n_years_model + n_years_proj x 2); which age,fleet of F to use for max F for msy/ypr calculations and projections
        spawn_seasons: n_stocks; which season spawning occurs for each stock
        spawn_regions: n_stocks; which region spawning occurs for each stock
        fleet_regions: n_fleets; which region each fleet is operating
        fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
             can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
           fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                   mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                    L: n_years_model x n_regions; "extra" mortality rate
               mature: n_stocks x n_years x n_ages; proportion mature
                  waa: (n_?) x n_years x n_ages_model. weight at age
      waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
       fracyr_seasons: n_seasons: length of intervals for each season
       avg_years_ind:
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_y = log_M_base.dim(2); 
  int n_ages = log_M_base.dim(3);
  array<Type> NAA_y(n_stocks, n_regions, n_ages);
  NAA_y.setZero();
  if(y==0) {
    NAA_y = get_NAA_1(N1_model,log_N1, NAA_where, log_M_base, FAA, which_F_age, 
      spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB, mu, L, mature, waa, waa_pointer_ssb, fracyr_seasons, 
      avg_years_ind, small_dim);
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
  matrix<Type> NAA_spawn_y(n_stocks,n_ages);
  NAA_spawn_y.setZero();
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++) {
    NAA_spawn_y(s,a) += NAA_y(s,r,a) * annual_SAA_spawn(s,y,a,r,spawn_regions(s)-1);
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
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions){
  /*
    provide "expected" recruitment (N(age 1)) for a given year
                  y: year (between 1 and n_years_model+n_years_proj)
      recruit_model: which recruitment model (1-4)
      mean_rec_pars: n_stocks x 2; recruitment parameters (defined in main code)
                SSB: n_years x n_stocks; of yearly SSB (uses y-1 for any S-R relationship)
                NAA: n_stocks x n_regions x n_years x n_ages; annual numbers at age by stock, region
           log_SR_a: yearly "a" parameters for SR function
           log_SR_b: yearly "b" parameters for SR function
         use_Ecov_R: 0/1 whether to use Ecov for recruitment
         Ecov_how_R: integer vector with an element that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: region where spawning and recruitment occur for each stock
  */
  int n_stocks = NAA.dim(0);
  vector<Type> pred_recruit(n_stocks);
  pred_recruit.setZero();
  for(int s = 0; s < n_stocks; s++){
    if(recruit_model(s) == 1) { // random walk
      pred_recruit(s) = NAA(s,spawn_regions(s)-1,y-1,0);
    } else {
      if(recruit_model(s) == 2) {// random about mean
        pred_recruit(s) = exp(mean_rec_pars(s,0));
        int nE = use_Ecov_R.rows(); 
        for(int i=0; i < nE; i++){
          if(use_Ecov_R(i,s) == 1) if(Ecov_how_R(i,s) == 1) pred_recruit(s) *= exp(Ecov_lm_R(s,y,i));
        }
      } else
      {
        if(recruit_model(s) == 3) // BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
        {
          pred_recruit(s) = exp(log_SR_a(y,s)) * SSB(y-1,s)/(1 + exp(log_SR_b(y,s))*SSB(y-1,s));
        } else // recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
        {
          pred_recruit(s) = exp(log_SR_a(y,s)) * SSB(y-1,s) * exp(-exp(log_SR_b(y,s)) * SSB(y-1,s));
        }
      }
    }
  }
  return pred_recruit;
}

template <class Type>
vector<Type> get_pred_recruit_y(int y, vector<int> recruit_model, matrix<Type> mean_rec_pars, vector<Type> SSB_y_minus_1, 
  array<Type> NAA_y_minus_1, matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions){
  /*
    provide "expected" recruitment (N(age 1)) for a given year
                  y: year (between 1 and n_years_model+n_years_proj)
      recruit_model: which recruitment model (1-4)
      mean_rec_pars: recruitment parameters (defined in main code)
                SSB_y_minus: n_stocks; SSB (uses y-1 for any S-R relationship) at previous year
                NAA_y_minus_1: n_stocks x n_regions x n_ages; numbers at age by stock, region at previous year
           log_SR_a: yearly "a" parameters for SR function
           log_SR_b: yearly "b" parameters for SR function
         use_Ecov_R: 0/1 whether to use Ecov for recruitment
         Ecov_how_R: integer vector with an element that tells how the Ecov is affecting recruitment
          Ecov_lm_R: (n_stocks, n_years_pop, n_Ecov); linear predictor for any environmental covariate effects on recruitment
      spawn_regions: region where spawning and recruitment occur for each stock
  */
  int n_stocks = NAA_y_minus_1.dim(0);
  vector<Type> pred_recruit(n_stocks);
  pred_recruit.setZero();
  for(int s = 0; s < n_stocks; s++){
    if(recruit_model(s) == 1) { // random walk
      pred_recruit(s) = NAA_y_minus_1(s,spawn_regions(s)-1,0);
    } else {
      if(recruit_model(s) == 2) {// random about mean
        pred_recruit(s) = exp(mean_rec_pars(s,0));
        int nE = use_Ecov_R.rows(); 
        for(int i=0; i < nE; i++){
          if(use_Ecov_R(i,s) == 1) if(Ecov_how_R(i,s) == 1) pred_recruit(s) *= exp(Ecov_lm_R(s,y,i));
        }
      } else
      {
        if(recruit_model(s) == 3) // BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
        {
          pred_recruit(s) = exp(log_SR_a(y,s)) * SSB_y_minus_1(s)/(1 + exp(log_SR_b(y,s))*SSB_y_minus_1(s));
        } else // recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
        {
          pred_recruit(s) = exp(log_SR_a(y,s)) * SSB_y_minus_1(s) * exp(-exp(log_SR_b(y,s)) * SSB_y_minus_1(s));
        }
      }
    }
  }
  return pred_recruit;
}

template <class Type>
array<Type> get_pred_N1(vector<int> N1_model, array<Type> log_N1, array<int> NAA_where, array<Type> N1_repars){
  /*
    provide the "expected" numbers at age in the first year. different from N1 only if N1 are random effects.
     N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
           N1: (n_stocks x n_regions x n_ages) parameters for estimated initial numbers at age
    NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
  */
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  array<Type> pred_N1(n_stocks, n_regions, n_ages);
  pred_N1.setZero();
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
    if(N1_model(s) < 2) {
      pred_N1(s,r,a) = exp(log_N1(s,r,a));
    } else{
      pred_N1(s,r,a) = exp(N1_repars(s,r,0)); //exp of mean of AR1 process
    } 
  }
  return(pred_N1);
}

template <class Type>
array<Type> get_pred_NAA_y(int y, vector<int> N1_model, array<Type> log_N1, array<Type> N1_repars, array<int> NAA_where, vector<int> recruit_model, 
  matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> Ps){

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
         use_Ecov_R: 0/1 whether to use Ecov for recruitment
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

  //array<Type> N1(n_stocks,n_regions,n_ages);
  if(y==0) { //initial NAA
    pred_NAA_y = get_pred_N1(N1_model, log_N1, NAA_where, N1_repars);
  } else {
  // Expected recruitment
    vector<Type> pred_recruit = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
      log_SR_b, use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions);
    for(int s = 0; s < n_stocks; s++){
      pred_NAA_y(s,spawn_regions(s)-1,0) = pred_recruit(s);
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
array<Type> get_pred_NAA_y(int y, vector<int> N1_model, array<Type> log_N1, array<Type> N1_repars, array<int> NAA_where, vector<int> recruit_model, 
  matrix<Type> mean_rec_pars, vector<Type> SSB_y_minus_1, array<Type> NAA_y_minus_1, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> Ps){

  /*
    provide "expected" numbers at age given NAA from previous time step
           N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
                 N1: n_stocks x n_regions x n_ages; array of parameters representing initial numbers at age.
          NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
      recruit_model: nstocks; which recruitment model (1-4)
      mean_rec_pars: n_stocks x 2; of any recruitment parameters (defined in main code)
                SSB_y_minus_1: n_stocks; of yearly SSB at previous year
                NAA_y_minus_1: array of numbers at age at previous year
           log_SR_a: yearly "a" parameters for SR function for each stock
           log_SR_b: yearly "b" parameters for SR function for each stock
         use_Ecov_R: 0/1 whether to use Ecov for recruitment
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

  //array<Type> N1(n_stocks,n_regions,n_ages);
  if(y==0) { //initial NAA
    pred_NAA_y = get_pred_N1(N1_model, log_N1, NAA_where, N1_repars);
  } else {
  // Expected recruitment
    vector<Type> pred_recruit = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB_y_minus_1, NAA_y_minus_1, log_SR_a, 
      log_SR_b, use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions);
    for(int s = 0; s < n_stocks; s++){
      pred_NAA_y(s,spawn_regions(s)-1,0) = pred_recruit(s);
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
array<Type> get_pred_NAA(int N1_model, array<Type> log_N1, array<Type> N1_repars, array<int> NAA_where, vector<int> recruit_model, 
  matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> annual_Ps){

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
         use_Ecov_R: 0/1 whether to use Ecov for recruitment
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
    array<Type> pred_NAA_y = get_pred_NAA_y(y, N1_model, log_N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB, NAA, 
      log_SR_a, log_SR_b, use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps);
    for(int a = 0; a < n_ages; a++) for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++){
      //pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
      pred_NAA(s,r,y,a) += pred_NAA_y(s,r,a);
    }
  }
  return(pred_NAA);
}

template <class Type>
array<Type> get_NAA(vector<int> NAA_re_model, vector<int> N1_model, array<Type> N1, array<Type> log_N1, array<Type> N1_repars, 
  array<Type> log_NAA, array<int> NAA_where, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb,
  vector<int> recruit_model, matrix<Type> mean_rec_pars, matrix<Type> log_SR_a, matrix<Type> log_SR_b, 
  matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> annual_Ps, array<Type> annual_SAA_spawn){
  /* 
    fill out numbers at age and "expected" numbers at age
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
             N1_model: 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations, 2: AR1 random effect
               N1: (n_stocks x n_regions x n_ages) numbers at age in the first year
               log_N1:
               N1_repars:
              log_NAA: (n_stocks x n_regions x nyears-1 x n_ages) parameters for ages after year 1
            NAA_where: n_stocks x n_regions x n_ages: 0/1 whether NAA exists in region at beginning of year. Also controls inclusion of any RE in nll.
               mature: n_stocks x n_years x n_ages; proportion mature
                  waa: (n_?) x n_years x n_ages_model. weight at age
      waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
      recruit_model:
      mean_rec_pars:
      log_SR_a:
      log_SR_b:
      use_Ecov_R:
      Ecov_how_R:
      Ecov_lm_R:
      spawn_regions:
      annual_Ps:
      annual_SAA_spawn:
  */
  int n_stocks = log_NAA.dim(0);
  int n_regions = log_NAA.dim(1);
  int n_y = log_NAA.dim(2)+1; 
  int n_ages = log_NAA.dim(3);
  array<Type> NAA(2,n_stocks, n_regions, n_y, n_ages); //NAA AND pred_NAA
  NAA.setZero();
  array<Type> NAA_last = N1;
  matrix<Type> NAA_spawn_last = get_NAA_spawn_y(0, NAA_last, annual_SAA_spawn, spawn_regions);
  vector<Type> SSB_last = get_SSB_y(0, NAA_spawn_last, waa, waa_pointer_ssb, mature);
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
    //NAA(s,r,0,a) = NAA1(s,r,a);
    NAA(0,s,r,0,a) = N1(s,r,a);
  }
  array<Type> pred_NAA_y = get_pred_NAA_y(0, N1_model, log_N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_last, NAA_last, 
    log_SR_a, log_SR_b, use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps);
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
    NAA(1,s,r,0,a) = pred_NAA_y(s,r,a);
  }

  for(int y = 1; y < n_y; y++){
    pred_NAA_y = get_pred_NAA_y(y, N1_model, log_N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_last, NAA_last, 
      log_SR_a, log_SR_b, use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps);
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) {
      NAA(1,s,r,y,a) = pred_NAA_y(s,r,a);
    }
    for(int s = 0; s < n_stocks; s++) {
      if(NAA_re_model(s) == 2){ //rec+1
        for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
          NAA(0,s,r,y,a) = exp(log_NAA(s,r,y-1,a));
        }
      }
      if(NAA_re_model(s) < 2) { //rec, Need to populate other ages with pred_NAA.
        NAA(0,s,spawn_regions(s)-1,y,0) = exp(log_NAA(s,spawn_regions(s)-1,y-1,0));
        for(int a = 1; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)){
          NAA(0,s,r,y,a) = pred_NAA_y(s,r,y-1,a);
        }
      }
      for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) NAA_last(s,r,a) = NAA(0,s,r,y,a);
      NAA_spawn_last = get_NAA_spawn_y(y, NAA_last,  annual_SAA_spawn, spawn_regions);
      SSB_last = get_SSB_y(y, NAA_spawn_last, waa, waa_pointer_ssb, mature);
    }
  }
  return(NAA);
}

template <class Type>
array<Type> extract_NAA(array<Type> all_NAA){
  int n_stocks = all_NAA.dim(1);
  int n_regions = all_NAA.dim(2);
  int n_y = all_NAA.dim(3); 
  int n_ages = all_NAA.dim(4);
  array<Type> NAA(n_stocks,n_regions,n_y,n_ages);
  NAA.setZero();
  for(int y = 1; y < n_y; y++) for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
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
  for(int y = 1; y < n_y; y++) for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    pred_NAA(s,r,y,a) = all_NAA(1,s,r,y,a);
  }
  return pred_NAA;
}

template <class Type>
array<Type> get_NAA_devs(array<Type> all_NAA){
  int n_stocks = all_NAA.dim(1);
  int n_regions = all_NAA.dim(2);
  int n_y = all_NAA.dim(3); 
  int n_ages = all_NAA.dim(4);
  array<Type> NAA_devs(n_stocks,n_regions,n_y,n_ages);
  NAA_devs.setZero();
  for(int y = 1; y < n_y; y++) for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++){
    NAA_devs(s,r,y,a) = log(all_NAA(0,s,r,y,a)) - log(all_NAA(1,s,r,y,a));
  }
  return NAA_devs;
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

  int n_y = Ecov_lm_R.dim(2);  
  int n_s = Ecov_lm_R.dim(0);
  int n_Ecov = Ecov_lm_R.dim(3);
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
  int n_y = Ecov_lm_R.dim(2);  
  int n_s = Ecov_lm_R.dim(0);
  int n_Ecov = Ecov_lm_R.dim(3);
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
array<Type> get_NAA_index(array<Type> NAA, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  matrix<Type> fracyr_indices, vector<int> index_seasons, vector<int> index_regions, array<Type> FAA, array<Type> log_M_base, 
  array<Type> mu, matrix<Type> L){
  /*
    produce the annual survival probabilities up to time of spawning for a given stock, age, season, year
                NAA: nstocks x nregions x nyears x nages; array of numbers at age 
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_seasons: n_seasons; length of intervals for each season
      fracyr_indices: n_indices; length of intervals for each index
      index_seasons: n_indices; which season the index occurs in
      index_regions: n_indices: which region the index is observing
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_indices = index_seasons.size();
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_years = log_M_base.dim(2);
  int n_ages = log_M_base.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim

  array<Type> NAA_index(n_stocks,n_indices,n_years,n_ages);
  NAA_index.setZero();
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) {
      for(int i = 0; i < n_indices; i++) {
        if(t == index_seasons(i)-1){ 
          //P(0,t) x P(t_i-t): PTM over interval from beginning of year to time of index within the season
          matrix<Type> P_index = P_y * get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_indices(y,i), 
            FAA, log_M_base, mu, L);
          for(int r = 0; r < n_regions; r++) NAA_index(s,i,y,a) += P_index(r,index_regions(i)-1) * NAA(s,r,y,a);
        }
      }
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      P_y = P_y * P_t;
    }
  }
  return(NAA_index);
}

template <class Type>
array<Type> get_NAA_catch(array<Type> NAA, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, 
  vector<Type> fracyr_seasons, array<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L){
  /*
    produce the numbers caught by stock, fleet, year, season, age up to time of spawning for a given stock, age, season, year
                NAA: nstocks x nregions x nyears x nages; array of numbers at age 
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons: 0/1 indicator whether fleet is operating in a given season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_seasons: n_seasons; length of intervals for each season
      fracyr_indices: n_indices; length of intervals for each index
      index_seasons: n_indices; which season the index occurs in
      index_regions: n_indices: which region the index is observing
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_years = log_M_base.dim(2);
  int n_ages = log_M_base.dim(3);
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
      matrix<Type> P_t = get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
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
matrix<Type> get_NAA_nll(array<Type> N1, vector<int> N1_model, array<Type> log_N1, array<Type> N1_repars, vector<int> NAA_re_model, array<Type> log_NAA, matrix<Type> log_NAA_sigma, matrix<Type> trans_NAA_rho, 
  array<int> NAA_where, vector<int> recruit_model, matrix<Type> mean_rec_pars, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R, 
  vector<int> spawn_regions, array<Type> annual_Ps, array<Type> annual_SAA_spawn, array<Type> waa, vector<int> waa_pointer_ssb, 
  array<Type> mature, int bias_correct_pe){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */
  //currently independent 2D (at most) AR1 processes by stock and region. not all ages may be available in all regions. number of stocks will typically be small
  //NAA_logsigma n_stocks x n_ages x n_regions
  //NAA_trans_rho n_stocks x 4 (rho_a, rho_y, rho_r, rho_s) 
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  array<Type> all_NAA = get_NAA(NAA_re_model, N1_model, N1, log_N1, N1_repars, log_NAA, NAA_where, 
    mature, waa, waa_pointer_ssb, recruit_model, mean_rec_pars, log_SR_a, log_SR_b, 
    use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions,  annual_Ps, annual_SAA_spawn); //log_NAA should be mapped accordingly to exclude NAA=0 e.g., recruitment by region.
  array<Type> NAA = extract_NAA(all_NAA);
  array<Type> pred_NAA = extract_pred_NAA(all_NAA);

  int n_stocks = NAA.dim(0);
  int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);
  int n_regions = NAA.dim(1);

  matrix<Type> marginal_sigma(n_ages,n_stocks);
  marginal_sigma.setZero();
  matrix<Type> nll_NAA(n_stocks,n_regions);
  nll_NAA.setZero();
  //can we do vector< vector< matrix<Type>>>?
  //NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0){
    // for NAA_re_model = 1, must make sure that rho_a = 0 and rho_y is set appropriately (cor = "iid" or "ar1_y") on R side
    Type NAA_rho_a = geninvlogit(trans_NAA_rho(s,0), Type(-1), Type(1), Type(2)); //2 is legacy
    Type NAA_rho_y = geninvlogit(trans_NAA_rho(s,1), Type(-1), Type(1), Type(2)); //2 is legacy
    marginal_sigma.col(s) = exp(vector<Type>(log_NAA_sigma.col(s))) * pow((1-pow(NAA_rho_y,2))*(1-pow(NAA_rho_a,2)),-0.5);
    if(NAA_re_model(s) == 1){ //"rec"
      vector<Type> NAA_devs_r_s(n_years-1);
      //for(int y = 1; y < n_years; y++) NAA_devs_r_s(y-1) = NAA_devs(s,spawn_regions(s)-1,y-1,0);
      for(int y = 1; y < n_years; y++) NAA_devs_r_s(y-1) = log(NAA(s,spawn_regions(s)-1,y,0)) - log(pred_NAA(s,spawn_regions(s)-1,y,0));
      if(bias_correct_pe) NAA_devs_r_s += 0.5*pow(marginal_sigma(0,s),2); //make sure this is ok when just recruitment is random.
      nll_NAA(s,spawn_regions(s)-1) += SCALE(AR1(NAA_rho_y), marginal_sigma(0,s))(NAA_devs_r_s);
    }
    if(NAA_re_model(s) == 2){ //"rec+1"
      for(int r = 0; r < n_regions; r++){
        int n_age_s_r = 0;
        for(int a = 0; a< n_ages; a++) if(NAA_where(s,r,a)) n_age_s_r++;
        if(n_age_s_r>0) {// has to be some fish of some age in this region
          array<Type> NAA_devs_s_r(n_years-1, n_age_s_r);
          NAA_devs_s_r.setZero();
          vector<Type> marginal_sigma_s_r(n_age_s_r);
          int k=0;
          for(int a = 0; a< n_ages; a++) if(NAA_where(s,r,a)) {
            marginal_sigma_s_r(k) = marginal_sigma(a,s);
            k++;
          }
          for(int y = 1; y < n_years; y++) {
            k=0;
            for(int a = 0; a< n_ages; a++) if(NAA_where(s,r,a)) {
              NAA_devs_s_r(y-1,k) = log(NAA(s,r,y,a)) - log(pred_NAA(s,r,y,a));
              if(bias_correct_pe) NAA_devs_s_r(y-1,k) += 0.5*pow(marginal_sigma_s_r(k),2);
              k++;
            }
          }
          nll_NAA(s,r) += SEPARABLE(VECSCALE(AR1(NAA_rho_a), marginal_sigma_s_r),AR1(NAA_rho_y))(NAA_devs_s_r);
        }
      }
    }
  }
  return nll_NAA;
}

template <class Type>
array<Type> simulate_NAA_devs(array<Type> NAA, vector<int> NAA_re_model, matrix<Type> log_NAA_sigma, matrix<Type> trans_NAA_rho, array<int> NAA_where, 
  vector<int> spawn_regions, int bias_correct_pe){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */

  //currently independent 2D (at most) AR1 processes by stock and region. not all ages may be available in all regions. number of stocks will typically be small
  //NAA_logsigma n_stocks x n_ages x n_regions
  //trans_NAA_rho n_stocks x 2 (rho_a, rho_y) 
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE

  int n_stocks = NAA.dim(0);
  int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);
  int n_regions = NAA.dim(1);
  array<Type> NAA_devs(n_stocks, n_regions, n_years-1, n_ages);
  matrix<Type> marginal_sigma(n_ages,n_stocks);
  marginal_sigma.setZero();

  //NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  for(int s = 0; s < n_stocks; s++) if(NAA_re_model(s)>0){
    // for NAA_re_model = 1, must make sure that rho_a = 0 and rho_y is set appropriately (cor = "iid" or "ar1_y") on R side
    Type NAA_rho_a = geninvlogit(trans_NAA_rho(s,0), Type(-1), Type(1), Type(2)); //2 is legacy
    Type NAA_rho_y = geninvlogit(trans_NAA_rho(s,1), Type(-1), Type(1), Type(2)); //2 is legacy
    marginal_sigma.col(s) = exp(vector<Type> (log_NAA_sigma.col(s))) * pow((1-pow(NAA_rho_y,2))*(1-pow(NAA_rho_a,2)),-0.5);
    if(NAA_re_model(s) == 1){ //"rec"
      vector<Type> NAA_devs_r_s(n_years-1);
      AR1(NAA_rho_y).simulate(NAA_devs_r_s); // sigma = 1, scale below
      NAA_devs_r_s = marginal_sigma(0,s) * NAA_devs_r_s;
      if(bias_correct_pe) NAA_devs_r_s -= 0.5*pow(marginal_sigma(0,s),2);
      for(int y = 0; y < n_years - 1; y++){
        NAA_devs(s,spawn_regions(s)-1,y,0) = NAA_devs_r_s(y);
      }
    }
    if(NAA_re_model(s) == 2){ //"rec+1"
      for(int r = 0; r < n_regions; r++){
        int n_age_s_r = 0;
        for(int a = 0; a< n_ages; a++) if(NAA_where(s,r,a)) n_age_s_r++;
        if(n_age_s_r>0) {// has to be some fish of some age in this region
          array<Type> NAA_devs_s_r(n_years-1, n_age_s_r);
          NAA_devs_s_r.setZero();
          vector<Type> marginal_sigma_s_r(n_age_s_r);
          int k=0;
          for(int a = 0; a< n_ages; a++) if(NAA_where(s,r,a)) marginal_sigma_s_r(k++) = marginal_sigma(a,s);
          SEPARABLE(VECSCALE(AR1(NAA_rho_a), marginal_sigma_s_r),AR1(NAA_rho_y)).simulate(NAA_devs_s_r); // scaled here
          for(int y = 1; y < n_years; y++) { 
            k=0;
            for(int a = 0; a< n_ages; a++) if(NAA_where(s,r,a)) {
              if(bias_correct_pe) NAA_devs_s_r(y-1,k) -= 0.5*pow(marginal_sigma_s_r(k),2);
              NAA_devs(s,r,y-1,a) = NAA_devs_s_r(y-1,k);
              k++;
            }
          }
        }
      }
    }
  }
  return NAA_devs;
}

// stopped here
template <class Type>
matrix<Type> simulate_log_NAA(vector<int> N1_model, array<Type> N1, array<Type> log_N1, array<Type> N1_repars, vector<int> NAA_re_model, array<Type> NAA_devs, array<Type> log_NAA,
  array<int> NAA_where, vector<int> recruit_model, matrix<Type> mean_rec_pars, 
  matrix<Type> log_SR_a, matrix<Type> log_SR_b, matrix<int> use_Ecov_R, matrix<int> Ecov_how_R, array<Type> Ecov_lm_R,
  vector<int> spawn_regions, array<Type> annual_Ps, array<Type> annual_SAA_spawn, array<Type> waa, vector<int> waa_pointer_ssb, 
  array<Type> mature){
  /*
            NAA_re_model: 0 SCAA, 1 "rec", 2 "rec+1"
  */
  //currently independent 2D (at most) AR1 processes by stock and region. not all ages may be available in all regions. number of stocks will typically be small
  //NAA_logsigma n_stocks x n_ages x n_regions
  //NAA_trans_rho n_stocks x 4 (rho_a, rho_y, rho_r, rho_s) 

  int n_stocks = log_NAA.dim(0);
  int n_regions = log_NAA.dim(1);
  int n_years = log_NAA.dim(2)+1;
  int n_ages = log_NAA.dim(3);
  array<Type> sim_log_NAA = log_NAA;
  sim_log_NAA.setZero();
  array<Type> NAA_y_minus_1 = N1;
  matrix<Type> NAA_spawn_y_minus_1 = get_NAA_spawn_y(0, NAA_y_minus_1,  annual_SAA_spawn, spawn_regions);
  vector<Type> SSB_y_minus_1 = get_SSB_y(0, NAA_spawn_y_minus_1, waa, waa_pointer_ssb, mature);
  
  for(int y = 1; y < n_years; y++) {
    //for y-1 = 0, it will populate with full initial NAA
    array<Type> pred_NAA_y = get_pred_NAA_y(y, N1_model, log_N1, N1_repars, NAA_where, recruit_model, mean_rec_pars, SSB_y_minus_1, NAA_y_minus_1, 
      log_SR_a, log_SR_b, use_Ecov_R, Ecov_how_R, Ecov_lm_R, spawn_regions, annual_Ps);
    NAA_y_minus_1.setZero();
    for(int s = 0; s < n_stocks; s++){
      if(NAA_re_model(s)==1){ //rec
        sim_log_NAA(s,spawn_regions(s)-1,y-1,0) += log(pred_NAA_y(s,spawn_regions(s)-1,0)) + NAA_devs(s,spawn_regions(s)-1,y-1,0);
        NAA_y_minus_1(s,spawn_regions(s)-1,0) = exp(sim_log_NAA(s,spawn_regions(s)-1,y-1,0));
        for(int a = 1; a < n_ages; a++) for(int r = 0; r < n_regions; r++) NAA_y_minus_1(s,r,a) = pred_NAA_y(s,r,a);
      }
      if(NAA_re_model(s)==2){ //rec+1
        for(int r = 0; r < n_regions; r++){
          for(int a = 0; a < n_ages; a++) if(NAA_where(s,r,a)) {
            sim_log_NAA(s,r,y-1,a) += log(pred_NAA_y(s,r,a)) + NAA_devs(s,r,y-1,a);
            NAA_y_minus_1(s,r,a) = exp(sim_log_NAA(s,r,y-1,a));
          }
        }
      }
    }
    NAA_spawn_y_minus_1 = get_NAA_spawn_y(y, NAA_y_minus_1, annual_SAA_spawn, spawn_regions);
    SSB_y_minus_1 = get_SSB_y(y, NAA_spawn_y_minus_1, waa, waa_pointer_ssb, mature);
  }

  return sim_log_NAA;
}
