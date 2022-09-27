template <class Type>
array<Type> get_eq_SAA(int y, vector<int> spawn_seasons, vector<int> fleet_regions, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_years x n_seasons x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: n_stocks; which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SAA.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,0,i,j) = S_ya(i,j);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I; //PTM for year and age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      } 
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,a,i,j) = S_ya(i,j);
    }
    //now plus group
    matrix<Type> P_ya = I; //PTM for year and age and up to time of spawning
    for(int t = 0; t < n_seasons; t++) {
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
    }
    matrix<Type> fundm(n_regions,n_regions);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1;
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,n_ages-1,i,j) = S_ya(i,j);
  }
  return SAA;
}

template <class Type>
array<Type> get_SPR(int y, vector<int> spawn_seasons, vector<int> fleet_regions, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
      fleet_regions: vector of indicators telling which region each fleet is operating
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
                FAA: fishing mortality: n_fleets x n_years x n_seasons x n_ages
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SPRAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SPRAA.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      matrix<Type> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        if(t == spawn_seasons(s)-1) {
          //P(0,t_spawn): PTM over entire year up to time of spawning
          P_spawn = P_ya * get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M_base, mu, L);
        }
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      } 
      // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
      //should be n_regions x n_regions
      matrix<Type> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,y,a) * waa(waa_pointer_ssb(s)-1,y,a);
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
        //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
        SPRAA(s,a,i,j) += SPR_ya(i,j); 
      }
      S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
    }
    //now plus group
    matrix<Type> P_ya = I, P_spawn = I; //PTM for year and age and up to time of spawning
    for(int t = 0; t < n_seasons; t++) {
      if(t == spawn_seasons(s)-1) {
        //P(0,t_spawn): PTM over entire year up to time of spawning
        P_spawn = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M_base, mu, L);
      }
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
    }
    matrix<Type> fundm(n_regions,n_regions);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1;
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // SSB/R at year and age = prob alive to up to age a-1 x prob spawn at age a x waa x mature
    //should be n_regions x n_regions
    matrix<Type> SPR_ya = S_ya * get_S(P_spawn, n_regions) * mature(s,y,n_ages-1) * waa(waa_pointer_ssb(s)-1,y,n_ages-1);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      //SSB per Recruit in each region (should only have positive values in spawn_regions(s)-1?)
      SPRAA(s,n_ages-1,i,j) += SPR_ya(i,j); 
    }
  }
  if(age_specific) {
    return SPRAA; 
  } else {
    array<Type> SPR(n_stocks,n_regions,n_regions);
    SPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      SPR(s,i,j) += SPRAA(s,a,i,j);
    }
    return SPR;
  }
}

template <class Type>
array<Type> get_SPR_0(int y, vector<int> spawn_seasons, array<int> can_move, 
  vector<int> mig_type, matrix<Type> fracyr_SSB, array<Type> log_M_base, array<Type> mu, matrix<Type> L, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_ssb, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate unfished equilibrium spawning biomass per recruit (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
       spawn_season: vector of indicators telling which season spawning occurs for each stock
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
         fracyr_SSB: n_years x n_stocks:  size of interval from beginning of season to time of spawning within that season
         log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
             mature: n_stocks x n_years x n_ages proportion mature at age
                waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_ssb: which waa matrix to use for ssb
     fracyr_seasons: n_seasons: length of intervals for each season
       age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_seasons = can_move.dim(2);
  int n_ages = can_move.dim(1);
  vector<int> fleet_regions(1);
  fleet_regions(0) = 1;
  array<Type> FAA(1,y+1,n_seasons,n_ages);
  FAA.setZero();

  array<Type> SPR_0 = get_SPR(y,spawn_seasons,fleet_regions,can_move,mig_type,fracyr_SSB,FAA,log_M_base,mu,L,mature,
    waa,waa_pointer_ssb, fracyr_seasons, age_specific, small_dim);
  
  return SPR_0;
}

template <class Type>
array<Type> get_YPR(int y, vector<int> fleet_regions, array<int> can_move, 
  vector<int> mig_type, array<Type> FAA, array<Type> log_M_base, array<Type> mu, matrix<Type> L, array<Type> waa, 
  vector<int> waa_pointer_fleet, vector<Type> fracyr_seasons, int age_specific, int small_dim){
  /* 
    calculate equilibrium yield per recruit (at age) by stock and region.
                    y: the model year for which to use SPR inputs
        fleet_regions: vector of indicators telling which region each fleet is operating
             can_move: n_stocks x ages x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
             mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                  FAA: fishing mortality: n_fleets x n_years x n_seasons x n_ages
           log_M_base: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                   mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                    L: n_years_model x n_regions. "extra" unobserved mortality
                  waa: (n_?) x n_years x n_ages_model. weight at age
    waa_pointer_fleet: which waa matrix to use for each fleet
       fracyr_seasons: n_seasons: length of intervals for each season
         age_specific: 0/1 telling whether to return SSB/R by age or not. dimensions of returned arrays are different. 
            small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M_base.dim(0);
  int n_regions = log_M_base.dim(1);
  int n_ages = log_M_base.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(2);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> YPRAA(n_stocks,n_ages,n_regions,n_fleets); //Yield/R at age in each fleet column, given recruited in region row
  YPRAA.setZero();
  matrix<Type> W(n_fleets,n_fleets);
  W.setZero();

  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < n_regions; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> S_ya = get_S(I, n_regions);
    for(int a = 0; a < n_ages-1; a++) {
      for(int i = 0; i < n_fleets; i++) W(i,i) = waa(waa_pointer_fleet(i)-1,y,a);
      matrix<Type> P_ya = I; //PTM for year and age
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P(a, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
      }
      // YPR at year and age = prob alive to up to age a-1 x prob caught at age a x waa
      matrix<Type> YPR_ya = S_ya * get_D(P_ya) * W; //should be n_regions x n_fleets
      S_ya = S_ya * get_S(P_ya, n_regions);  //accumulate for next age
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,a,i,j) = YPR_ya(i,j);
    }
    //now plus group
    for(int i = 0; i < n_fleets; i++) W(i,i) = waa(waa_pointer_fleet(i)-1,y,n_ages-1);
    matrix<Type> P_ya = I; //PTM for year and age
    for(int t = 0; t < n_seasons; t++) {
      //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
      P_ya = P_ya * get_P(n_ages-1, y, s, t, fleet_regions, can_move, mig_type, fracyr_seasons(t), FAA, log_M_base, mu, L);
    }
    matrix<Type> fundm(n_regions,n_regions);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
      fundm(i,j) = -P_ya(i,j);
      if(i==j) fundm(i,j) += 1;
    }
    if(small_dim) fundm = fundm.inverse(); else fundm = matinv(fundm);
    //for plus group S_ya = S_y,a-1 x (I - S_y,+)^-1
    S_ya = S_ya * fundm;
    // YPR at year and age = (prob alive to up to age a-1 + prob alive at older ages) x prob caught at these ages x waa
    matrix<Type> YPR_ya = S_ya * get_D(P_ya) * W; //should be n_regions x n_fleets
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) YPRAA(s,n_ages-1,i,j) = YPR_ya(i,j);
  }
  if(age_specific) {
    return YPRAA; 
  } else {
    array<Type> YPR(n_stocks,n_regions,n_fleets);
    YPR.setZero();
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) {
      YPR(s,i,j) += YPRAA(s,a,i,j);
    }
    return YPR;
  }
}

/* calculate single SSB/R at F for spatial model across stocks and regions */

template<class Type>
struct spr_F_spatial {
  // Data and parameter objects for calculation: 
  vector<int> spawn_seasons;
  vector<int> fleet_regions; 
  array<int> can_move; 
  vector<int> mig_type;
  matrix<Type> fracyr_SSB;
  array<Type> selectivityAA;
  array<Type> log_M_base;
  array<Type> mu; 
  matrix<Type> L;
  array<Type> mature;
  array<Type> waa; 
  vector<int> waa_pointer_ssb;
  vector<Type> fracyr_seasons;
  int age_specific; 
  int small_dim;
  int y;

  // Constructor 
  spr_F_spatial(
  vector<int> spawn_seasons_,
  vector<int> fleet_regions_, 
  array<int> can_move_, 
  vector<int> mig_type_,
  matrix<Type> fracyr_SSB_,
  array<Type> selectivityAA_,
  array<Type> log_M_base_,
  array<Type> mu_, 
  matrix<Type> L_,
  array<Type> mature_,
  array<Type> waa_, 
  vector<int> waa_pointer_ssb_,
  vector<Type> fracyr_seasons_,
  int age_specific_, 
  int small_dim_,
  int y_) :
    spawn_seasons(spawn_seasons_), 
    fleet_regions(fleet_regions_),
    can_move(can_move_),
    mig_type(mig_type_),
    fracyr_SSB(fracyr_SSB_),
    selectivityAA(selectivityAA_),
    log_M_base(log_M_base_),
    mu(mu_),
    L(L_),
    mature(mature_),
    waa(waa_),
    waa_pointer_ssb(waa_pointer_ssb_),
    fracyr_seasons(fracyr_seasons_),
    age_specific(age_specific_),
    small_dim(small_dim_),
    y(y_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_stocks = log_M_base.dim(0);
    int n_regions = log_M_base.dim(1);
    int n_ages = log_M_base.dim(3);
    int n_fleets = fleet_regions.size();
    int n_seasons = can_move.dim(2);
    array<Type> FAA(n_fleets,y+1,n_seasons,n_ages);
    for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) {
      FAA(f,y,t,a) = selectivityAA(f,y,t,a) * exp(Type(log_F));
    }
    array<T> SPR_sr = get_SPR(y, spawn_seasons, fleet_regions, can_move, mig_type, fracyr_SSB, FAA, log_M_base, mu, L, mature, waa, 
      waa_pointer_ssb, fracyr_seasons, 0, small_dim).template cast<T>();

    T SPR = SPR_sr.sum();
    return SPR;
  }
};

