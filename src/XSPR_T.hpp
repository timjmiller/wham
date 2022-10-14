template <class T>
matrix<T> get_P(int age, int year, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, T time, array<T> FAA, array<T>log_M, 
  array<T> mu, matrix<T> L) {
  /*
    produce the probability transition matrix for a given stock,, age, season, year
                age: which age
               year: which year
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_regions = log_M.dim(1);
  int n_fleets = fleet_regions.size();
  int dim = n_regions+n_fleets+1;
  matrix<T> P(dim,dim);
  P.setZero(); //zero it out.
  vector<T> F(n_fleets), M(n_regions), Z(n_regions);
  F.setZero();
  M.setZero();
  Z.setZero();
  for(int r = 0; r < n_regions; r++) {
    M(r) = exp(log_M(stock,r,year,age)) + L(year,r); //add in "extra" mortality
    Z(r) = M(r);
  }
  for(int f = 0; f < n_fleets; f++) {
    if(fleet_seasons(f,season)) F(f) = FAA(year,f,age);
    Z(fleet_regions(f)-1) += F(f);
  }
  if(n_regions == 1){ //usual Baranov
    if(time < 1e-15) {
      P(0,0) = 1.0;
    } else {
      P(0,0) = exp(-Z(0)* time);
      for(int f = 0; f < n_fleets; f++) {
        P(0,1+f) = F(f) * (1 - exp(-Z(0) * time))/Z(0);
      }
      P(0,1+n_fleets) = M(0) * (1 - exp(-Z(0) * time))/Z(0);
    }
    for(int i = 1; i < dim; i++) P(i,i) = 1.0;
  } else{ //n_regions > 1
    matrix<int> can_move_now(n_regions,n_regions);
    can_move_now.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) can_move_now(i,j) = can_move(stock,age,season,i,j);
    //if(n_mu(stock,year,season,age)>0) //migration is happening
    if(can_move_now.sum()>0) { //migration is happening
      if(mig_type(stock) == 0) { //migration is instantaneous after survival and mortality, so P is easy.
        if(time < 1e-15) {
        //prob of survival is 1 when interval is zero
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) P(i,j) =  mu(stock,age,season,year,i,j);
        } else {
          for(int f = 0; f < n_fleets; f++) P(fleet_regions(f)-1,n_regions + f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
          for(int i = 0; i < n_regions; i++) {
            //probs of survival and mortality and then moving between regions
            for(int j = 0; j < n_regions; j++) 
            {
              //survival only a function of region starting in, then modify survival in in each region by prob of moving.
              P(i,j) = exp(-Z(i) * time) * mu(stock,age,season,year,i,j); 
            }
            //caught only in region starting in since migration happens after survival 
            
            //all other dead
            P(i,dim-1) = M(i) * (1.0 - exp(-Z(i) * time))/Z(i);
          }
        }
        //fish and nat mort state: if dead, must stay dead
        for(int i = n_regions; i < dim; i++) P(i,i) = 1.0; 
      }
      if(mig_type(stock) == 1) {//migration occurs continuously during interval, so P is not so easy.
        //prob of survival is 1 when interval is zero
        if(time < 1e-15) {
          for(int i = 0; i < n_regions; i++) P(i,i) = 1.0;
        } else {
          matrix<T> A(dim,dim);
          A.setZero(); //zero it out.
          for(int i = 0; i < n_regions; i++) A(i,dim-1) += M(i); //other dead
          for(int f = 0; f < n_fleets; f++) A(fleet_regions(f)-1,n_regions + f) += F(f);
          for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
            //need to not put in mu diagonal because it has -sum of mig rates in there already.
            if(i != j) A(i,j) += mu(stock,age,season,year,i,j); // transition intensities
          }
          for(int i = 0; i < n_regions; i++) A(i,i) = -(A.row(i)).sum(); //hazard
          A = A * time;
          P = expm(A);
        }
      }
    }
    else {//no migration during this interval, so P is easy.
      //prob of survival is 1 when interval is zero
      if(time < 1e-15) {
        for(int i = 0; i < dim; i++) P(i,i) = 1.0; 
      } else {
        for(int i = 0; i < n_regions; i++) 
        {
          P(i,i) = exp(-Z(i) * time);
          //other dead
          P(i,dim-1) = M(i) * (1.0 - exp(-Z(i) * time))/Z(i);
        }
        for(int f = 0; f < n_fleets; f++) 
        {
          P(fleet_regions(f)-1,n_regions+f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
        } 
      }
    }
  } //end n_regions > 1
  return(P);
}

template <class T>
matrix<T> get_P(int age, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, T time, matrix<T> FAA, array<T>log_M, 
  array<T> mu, vector<T> L) {
  /*
    produce the probability transition matrix for a given stock, age, season
                age: which age
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions; movement rates
                  L: n_regions; "extra" mortality rate
  */
  array<T> FAA_d(FAA.rows(),1, FAA.cols());
  FAA_d.setZero();
  for(int f = 0; f < FAA.rows(); f++) for(int a = 0; a < FAA.cols(); a++) FAA_d(f,0,a) = FAA(f,a);
  array<T> log_M_d(log_M.dim(0),log_M.dim(1),1,log_M.dim(2));
  log_M_d.setZero();
  for(int s = 0; s < log_M.dim(0); s++) for(int a = 0; a < log_M.dim(2); a++) for(int r = 0; r < log_M.dim(1); r++) {
    log_M_d(s,r,0,a) = log_M(s,r,a);
  }
  array<T> mu_d(mu.dim(0), mu.dim(1), mu.dim(2), 1, mu.dim(3),mu.dim(4));
  mu_d.setZero();
  for(int s = 0; s < mu.dim(0); s++) for(int a = 0; a < mu.dim(1); a++) for(int t = 0; t < mu.dim(2); t++) {
    for(int r = 0; r < mu.dim(3); r++) for(int rr = 0; rr < mu.dim(4); rr++){
      mu_d(s,a,t,0,r,rr) = mu_d(s,a,t,r,rr);
    }
  }
  matrix<T> L_d(1,L.size());
  L_d.setZero();
  for(int r = 0; r < mu.dim(3); r++) L_d(0,r) = L(r);
  matrix<T> P = get_P(age, 0, stock, season, fleet_regions, fleet_seasons, can_move, mig_type, time, FAA_d, log_M_d, mu_d, L_d);
  return(P);
}

template <class T>
matrix<T> get_S(matrix<T> P, int n_regions){
  /*
    extract the submatrix from a PTM that contains the proportions surviving in each region
              P: the probablity transition matrix
      n_regions: the number of regions
  */
  matrix<T> S(n_regions,n_regions);
  S.setZero();
  for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) S(i,j) = P(i,j);
  return(S);
}


template <class T>
array<T> get_SPR(vector<int> spawn_seasons, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, vector<T> fracyr_SSB, matrix<T> FAA, array<T> log_M, array<T> mu, vector<T> L, 
  matrix<T> mature, matrix<T> waa_ssb, vector<T> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
    NB: This version of get_SPR does not allow year indices. It is assumed that the yearly (or averaged) values are provided.
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions array of movement matrices
                  L: n_regions. "extra" unobserved mortality
             mature: n_stocks x n_ages proportion mature at age
            waa_ssb: n_stocks x n_ages. weight at age
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_ages = log_M.dim(3);
  int n_fleets = FAA.rows();
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<T> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();

  matrix<T> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<T> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<T> P_ya = I, P_spawn = I; //PTM for age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        if(t == spawn_seasons(s)-1) {
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M, mu, L);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      } 
      // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
      //should be n_regions x n_regions
      matrix<T> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,a) * waa_ssb(s,a);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
        SPRAA(s,a,i,j) += SPR_ya(i,j); 
      }
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
    }
    //now plus group
    matrix<T> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
    for(int t = 0; t < n_seasons; t++) {
      if(t == spawn_seasons(s)-1) {
        //P(0,t_spawn): PTM over entire year up to time of spawning
        P_spawn = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(s), FAA, log_M, mu, L);
      }
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    }
    matrix<T> fundm(n_regions,n_regions);
    fundm.setZero();
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1;
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
    //should be n_regions x n_regions
    matrix<T> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,n_ages-1) * waa_ssb(s,n_ages-1);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
      SPRAA(s,n_ages-1,i,j) += SPR_ya(i,j); 
    }
  }
  if(age_specific) {
    return SPRAA; 
  } else {//sum over ages
    array<T> SPR(n_stocks,n_regions,n_regions);
    SPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      SPR(s,i,j) += SPRAA(s,a,i,j);
    }
    return SPR;
  }
}

/* calculate single SSB/R at F for spatial model across stocks and regions */

template<class Type>
struct spr_F_spatial {
  // Data and parameter objects for calculation: 
  vector<int> spawn_seasons;
  vector<int> spawn_regions;
  vector<int> fleet_regions; 
  matrix<int> fleet_seasons; 
  array<int> can_move; 
  vector<int> mig_type;
  vector<Type> fracyr_SSB;
  matrix<Type> selectivity;
  array<Type> log_M;
  array<Type> mu; 
  vector<Type> L;
  matrix<Type> mature;
  matrix<Type> waa_ssb; 
  vector<Type> fracyr_seasons;
  vector<Type> SPR_weights; //how to weight stock-specific SSB/R for aggregate SSB/R.
  int age_specific; 
  int small_dim;

  // Constructor 
  spr_F_spatial(
  vector<int> spawn_seasons_,
  vector<int> spawn_regions_,
  vector<int> fleet_regions_, 
  matrix<int> fleet_seasons_, 
  array<int> can_move_, 
  vector<int> mig_type_,
  vector<Type> fracyr_SSB_,
  matrix<Type> selectivity_,
  array<Type> log_M_,
  array<Type> mu_, 
  vector<Type> L_,
  matrix<Type> mature_,
  matrix<Type> waa_ssb_, 
  vector<Type> fracyr_seasons_,
  vector<Type> SPR_weights_,
  int age_specific_, 
  int small_dim_) :
    spawn_seasons(spawn_seasons_), 
    spawn_regions(spawn_regions_), 
    fleet_regions(fleet_regions_),
    fleet_seasons(fleet_seasons_),
    can_move(can_move_),
    mig_type(mig_type_),
    fracyr_SSB(fracyr_SSB_),
    selectivity(selectivity_),
    log_M(log_M_),
    mu(mu_),
    L(L_),
    mature(mature_),
    waa_ssb(waa_ssb_),
    fracyr_seasons(fracyr_seasons_),
    SPR_weights(SPR_weights_),
    age_specific(age_specific_),
    small_dim(small_dim_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_stocks = log_M.dim(0);
    int n_regions = log_M.dim(1);
    int n_ages = log_M.dim(3);
    int n_fleets = fleet_regions.size();
    int n_seasons = can_move.dim(2);
    matrix<T> FAA(n_fleets,n_ages);
    FAA.setZero();
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) {
      FAA(f,a) = T(selectivity(f,a)) * exp(log_F(0));
    }
    vector<T> fracyrssbT = fracyr_SSB.template cast<T>();
    array<T> logMbaseT(log_M.dim(0),log_M.dim(1),log_M.dim(2));
    for(int i = 0; i < log_M.dim(0); i++) for(int j = 0; j < log_M.dim(1); j++) for(int k = 0; k < log_M.dim(2); k++){
      logMbaseT(i,j,k) = T(log_M(i,j,k));
    }
    array<T> muT(mu.dim(0),mu.dim(1),mu.dim(2),mu.dim(3),mu.dim(4));
    for(int i = 0; i < mu.dim(0); i++) for(int j = 0; j < mu.dim(1); j++) for(int k = 0; k < mu.dim(2); k++){
      for(int l = 0; l < mu.dim(0); l++) for(int m = 0; m < mu.dim(4); m++) {
        muT(i,j,k,l,m) = T(mu(i,j,k,l,m));
      }
    }
    vector<T> LT = L.template cast<T>();
    matrix<T> matT = mature.template cast<T>();
    matrix<T> waassbT = waa_ssb.template cast<T>();
    vector<T> fracyrseasonT = fracyr_seasons.template cast<T>();
    //stock-specific SSB/R
    array<T> SPR_sr = get_SPR(spawn_seasons, fleet_regions, fleet_seasons, can_move, mig_type, 
      fracyrssbT, 
      FAA, 
      logMbaseT, 
      muT, 
      LT, 
      matT, 
      waassbT, 
      fracyrseasonT, 
      0, small_dim);

    //weighted-average of SSB/R returned
    T SPR = 0;
    for(int s = 0; s < n_stocks; s++) SPR += T(SPR_weights(s)) * SPR_sr(s,spawn_regions(s)-1,spawn_regions(s)-1); 
    
    return SPR;
  }
};
