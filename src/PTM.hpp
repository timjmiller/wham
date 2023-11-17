//NOTE get_P_t_base here is defined as class T instead of Type, but is currently used interchangeably.
// Not sure if this affects expected model performance.
template <class T>
matrix<T> get_P_t_base(vector<int> fleet_regions, matrix<int> can_move, int mig_type, T time, vector<T> F, vector<T> M, 
  matrix<T> mu, vector<T> L, int trace = 0) {
  /*
    produce the probability transition matrix over a time interval
      fleet_regions: n_fleets; which region each fleet is operating
           can_move: n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                  F: fishing mortality: n_fleets
                  M: n_regions
                 mu: n_regions x n_regions; movement rates
                  L: n_regions; "extra" mortality rate
              trace: 0/1 whether to print stuff to screen
  */
  int n_regions = L.size();
  int n_fleets = fleet_regions.size();
  int dim = n_regions+n_fleets+1;
  if(trace) see("inside get_P_t_base");
  if(trace) see(F);
  if(trace) see(M);
  matrix<T> P(dim,dim);
  P.setZero(); //zero it out.
  vector<T> Z(n_regions);
  Z.setZero();
  for(int r = 0; r < n_regions; r++) {
    //M(r) = exp(log_M(stock,r,year,age)) + L(year,r); //add in "extra" mortality
    Z(r) = M(r) + L(r);
  }
  if(trace) see(Z);
  for(int f = 0; f < n_fleets; f++) {
    Z(fleet_regions(f)-1) += F(f);
  }
  if(trace) see(Z);
  if(trace) see(time);
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
    if(trace) see(P);
  } else{ //n_regions > 1
    //matrix<int> can_move_now(n_regions,n_regions);
    //can_move_now.setZero();
    //for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) can_move_now(i,j) = can_move(stock,season,i,j);
    if(can_move.sum()>0) { //migration is happening
      if(mig_type == 0) { //migration is instantaneous after survival and mortality, so P is easy.
        if(time < 1e-15) {
        //prob of survival is 1 when interval is zero: P = I * P_move
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) P(i,j) =  mu(i,j);
        } else {
          for(int f = 0; f < n_fleets; f++) P(fleet_regions(f)-1,n_regions + f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
          for(int i = 0; i < n_regions; i++) {
            //probs of survival and mortality and then moving between regions
            for(int j = 0; j < n_regions; j++) 
            {
              //survival only a function of region starting in, then modify survival in in each region by prob of moving.
              P(i,j) = exp(-Z(i) * time) * mu(i,j); 
            }
            //caught only in region starting in since migration happens after survival 
            
            //all other dead
            P(i,dim-1) = (M(i)+L(i)) * (1.0 - exp(-Z(i) * time))/Z(i);
          }
        }
        //fish and nat mort state: if dead, must stay dead
        for(int i = n_regions; i < dim; i++) P(i,i) = 1.0; 
      }
      if(mig_type == 1) {//migration occurs continuously during interval, so P is not so easy.
        //prob of survival and staying is 1 when interval is zero
        if(time < 1e-15) {
          for(int i = 0; i < n_regions; i++) P(i,i) = 1.0;
        } else {
          matrix<T> A(dim,dim);
          A.setZero(); //zero it out.
          for(int i = 0; i < n_regions; i++) A(i,dim-1) += M(i) + L(i); //other dead
          for(int f = 0; f < n_fleets; f++) A(fleet_regions(f)-1,n_regions + f) += F(f);
          for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
            //need to not put in mu diagonal because it has -sum of mig rates in there already.
            if(i != j) A(i,j) += mu(i,j); // transition intensities
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
          P(i,dim-1) = (M(i) + L(i)) * (1.0 - exp(-Z(i) * time))/Z(i);
        }
        for(int f = 0; f < n_fleets; f++) 
        {
          P(fleet_regions(f)-1,n_regions+f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
        } 
        for(int i = n_regions; i < dim; i++) P(i,i) = 1.0; 
      }
    }
  } //end n_regions > 1
  if(trace) see("end inside get_P_t_base");

  return P;
}

//NOTE get_P_t here is defined as class T instead of Type and is not distiguishabled when used.
//Not sure if this affects expected model performance.
template <class T>
matrix<T> get_P_t(int age, int year, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, T time, array<T> FAA, array<T>log_M, 
  array<T> mu, matrix<T> L, int trace = 0) {
  /*
    produce the probability transition matrix for a given stock, age, season, year
                age: which age
               year: which year
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
              trace: 0/1 whether to print stuff to screen
  */
  int n_regions = log_M.dim(1);
  //int n_fleets = fleet_regions.size();
  if(trace) see("inside get_P_t");
  vector<T> M(n_regions);
  vector<T> F = get_F_t(vector<int> (fleet_seasons.col(season)), age, year, FAA);
  if(trace) see(F);
  matrix<T> mu_stya(n_regions,n_regions);
  M.setZero();
  mu_stya.setZero();
  for(int r = 0; r < n_regions; r++) M(r) = exp(log_M(stock,r,year,age));
  if(trace) see(M);
  matrix<int> can_move_sta(n_regions, n_regions);
  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
    can_move_sta(r,rr) = can_move(stock,season,r,rr);
    mu_stya(r,rr) = mu(stock,age,season,year,r,rr);
  }
  if(trace) see(can_move_sta);
  if(trace) see(mu_stya);
  if(trace) see(time);
  matrix<T> P = get_P_t_base(fleet_regions, can_move_sta, mig_type(stock), time, F, M, mu_stya, vector<T> (L.row(year)), trace);
  if(trace) see(P);
  if(trace) see("end inside get_P_t");
  return P;
}

template <class T>
matrix<T> get_P_t(int age, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, T time, matrix<T> FAA, array<T>log_M, 
  array<T> mu, vector<T> L, int trace = 0) {
  /*
    produce the probability transition matrix for a given stock, age, season FROM YEAR-SPECIFIC PARAMETERS
                age: which age
              stock: which stock
             season: which season
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
               time: time interval between beginning and end of the interval (units = years)
                FAA: fishing mortality: n_fleets x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_regions x n_regions; movement rates
                  L: n_regions; "extra" mortality rate
              trace: 0/1 whether to print stuff to screen
  */
  int n_regions = log_M.dim(1);
  //int n_fleets = fleet_regions.size();
  if(trace) see("inside get_P_t from year-specific parameters");
  vector<T> M(n_regions);
  vector<T> F = get_F_t(vector<int> (fleet_seasons.col(season)), age, FAA);
  if(trace) see(F);
  matrix<T> mu_stya(n_regions,n_regions);
  M.setZero();
  mu_stya.setZero();
  for(int r = 0; r < n_regions; r++) M(r) = exp(log_M(stock,r,age));
  if(trace) see(M);
  matrix<int> can_move_sta(n_regions, n_regions);
  for(int r = 0; r < n_regions; r++) for(int rr = 0; rr < n_regions; rr++) {
    can_move_sta(r,rr) = can_move(stock,season,r,rr);
    mu_stya(r,rr) = mu(stock,age,season,r,rr);
  }
  if(trace) see(can_move_sta);
  if(trace) see(mu_stya);
  if(trace) see(time);
  matrix<T> P = get_P_t_base(fleet_regions, can_move_sta, mig_type(stock), time, F, M, mu_stya, L, trace);
  if(trace) see(P);
  if(trace) see("end inside get_P_t");
  return P;
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
matrix<T> get_D(matrix<T> P, int n_regions, int n_fleets){
  /*
    extract the submatrix from a PTM that contains the proportions captured by each fleet in each region
              P: the probablity transition matrix
      n_regions: the number of regions
       n_fleets: the number of fleets
  */
  matrix<T> D(n_regions,n_fleets);
  D.setZero();
  for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) D(i,j) = P(i,j+n_regions);
  return(D);
}

template <class Type>
array<Type> get_annual_Ps(int n_years_model, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the annual probability transition matrix for a given stock, age, season, year
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_seasons: n_seasons; length of intervals for each season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions; movement rates
                  L: n_years x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_years = log_M.dim(2);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim
  //get probability transition matrices for yearly survival, movement, capture...
  //also get annual NAA at spawning and NAA for each index along the way.
  array<Type> annual_Ps(n_stocks,n_years,n_ages,P_dim,P_dim);
  annual_Ps.setZero();  
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();  
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) {
      //get numbers at age a for stock s in each region at time of spawning
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
    for(int i = 0; i < P_dim; i++) for(int j = 0; j < P_dim; j++) annual_Ps(s,y,a,i,j) = P_y(i,j);
  }
  return annual_Ps;
}

template <class Type>
array<Type> update_annual_Ps(int y, array<Type> annual_Ps, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the annual probability transition matrix for a given stock, age, season, year
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_seasons: n_seasons; length of intervals for each season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions; movement rates
                  L: n_years x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  //int n_years = log_M.dim(2);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim
  //get probability transition matrices for yearly survival, movement, capture...
  //also get annual NAA at spawning and NAA for each index along the way.
  array<Type> updated_annual_Ps = annual_Ps;
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();  
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) {
      //get numbers at age a for stock s in each region at time of spawning
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
    for(int i = 0; i < P_dim; i++) for(int j = 0; j < P_dim; j++) updated_annual_Ps(s,y,a,i,j) = P_y(i,j);
  }
  return updated_annual_Ps;
}

template <class Type>
array<Type> get_annual_SAA_spawn(int n_years_model, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  matrix<Type> fracyr_SSB, vector<int> spawn_seasons, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the annual survival probabilities up to time of spawning for a given stock, age, season, year
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_sesons: n_seasons; length of intervals for each season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_years = log_M.dim(2);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim

  array<Type> annual_SAA_SSB(n_stocks,n_years,n_ages,n_regions,n_regions); //just survival categories
  annual_SAA_SSB.setZero();  
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();  
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < spawn_seasons(s)-1; t++) { // only need to go up to the season prior to when spawning occurs
      //get numbers at age a for stock s in each region at time of spawning
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
    //P(0,t) x P(t_s-t): PTM over interval from to time of spawning within the season
    matrix<Type> P_SSB = P_y * get_P_t(a, y, s, spawn_seasons(s)-1, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M, mu, L);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) annual_SAA_SSB(s,y,a,i,j) = P_SSB(i,j);
  }
  return annual_SAA_SSB;
}

template <class Type>
array<Type> update_annual_SAA_spawn(int y, array<Type> annual_SAA_spawn, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  matrix<Type> fracyr_SSB, vector<int> spawn_seasons, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the annual survival probabilities up to time of spawning for a given stock, age, season, year
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_sesons: n_seasons; length of intervals for each season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years_pop x n_regions x n_regions; movement rates
                  L: n_years_model x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  //int n_years = log_M.dim(2);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim

  array<Type> updated_annual_SAA_spawn = annual_SAA_spawn;
  matrix<Type> I_mat(P_dim,P_dim);
  I_mat.setZero();  
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;
  for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < spawn_seasons(s)-1; t++) { // only need to go up to the season prior to when spawning occurs
      //get numbers at age a for stock s in each region at time of spawning
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
    //P(0,t) x P(t_s-t): PTM over interval from to time of spawning within the season
    matrix<Type> P_SSB = P_y * get_P_t(a, y, s, spawn_seasons(s)-1, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M, mu, L);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) updated_annual_SAA_spawn(s,y,a,i,j) = P_SSB(i,j);
  }
  return updated_annual_SAA_spawn;
}

template <class Type>
array<Type> get_seasonal_Ps_y(int y, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, 
  vector<Type> fracyr_seasons, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the probability transition matrices for each stock, season, age for year y
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous
      fracyr_sesons: n_seasons; length of intervals for each season
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions; movement rates
                  L: n_years x n_regions; "extra" mortality rate
  */
  int n_fleets = FAA.dim(0);
  int n_seasons = fleet_seasons.cols();
  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_ages = log_M.dim(3);
  int P_dim = n_regions + n_fleets + 1; // probablity transition matrix is P_dim x P_dim
  array<Type> P_seasonal_y(n_stocks,n_seasons,n_ages,P_dim,P_dim);
  P_seasonal_y.setZero();
  for(int s = 0; s < n_stocks; s++) for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) 
  {
    matrix<Type> P_t = get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    for(int d = 0; d < P_dim; d++) for(int dd = 0; dd < P_dim; dd++) P_seasonal_y(s,t,a,d,dd) = P_t(d,dd);
  }
  return P_seasonal_y;
}


template <class Type>
array<Type> get_eq_SAA(int y, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, 
  vector<int> mig_type, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L, 
  vector<Type> fracyr_seasons, int small_dim){
  /* 
    calculate equilibrium survival (at age) by stock and region. If movement is set up approriately 
    all fish can be made to return to a single spawning region for each stock.
                  y: the model year for which to use SPR inputs
      fleet_regions: vector of indicators telling which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x n_seasons x n_regions x n_regions: 0/1 determining whether movement can occur from one region to another
           mig_type: n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
                FAA: fishing mortality: n_fleets x n_years x n_ages
         log_M: log M (density-independent components): n_stocks x n_regions x ny x n_ages
                 mu: n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions array of movement matrices
                  L: n_years_model x n_regions. "extra" unobserved mortality
     fracyr_seasons: n_seasons: length of intervals for each season
          small_dim: 0/1 telling whether the n_regions is "small." Different methods of inverting matrices.
  */

  int n_stocks = log_M.dim(0);
  int n_regions = log_M.dim(1);
  int n_ages = log_M.dim(3);
  int n_fleets = FAA.dim(0);
  int n_seasons = can_move.dim(1);
  int P_dim = n_regions + n_fleets + 1;
  array<Type> SAA(n_stocks,n_ages,n_regions,n_regions); //SSB/R at age in each region column, given recruited in region row
  SAA.setZero();
  matrix<Type> I(P_dim,P_dim);
  I.setZero();
  for(int i = 0; i < P_dim; i++) I(i,i) = 1.0;

  for(int s = 0; s < n_stocks; s++) {
    matrix<Type> P_ya = I; //PTM for year and age and up to time of spawning
    matrix<Type> S_ya = get_S(P_ya, n_regions);
    for(int a = 0; a < n_ages; a++) {
      matrix<Type> P_ya = I; //PTM for year and age and up to time of spawning
      for(int t = 0; t < n_seasons; t++) {
        //update PTM to end of season t P(0,s) * P(s,t) = P(0,t)
        P_ya = P_ya * get_P_t(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      }
      if(a == n_ages-1){
        //plus group
        matrix<Type> fundm(n_regions,n_regions);
        fundm.setZero();
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
          fundm(i,j) = -P_ya(i,j);
          if(i==j) fundm(i,j) += 1;
        }
        if(small_dim) fundm = fundm.inverse(); else fundm = atomic::matinv(fundm);
        //for plus group S_ya = cum(S_y,a-1) x (I - S_y,+)^-1
        S_ya = S_ya * fundm;
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,a,i,j) = S_ya(i,j);
      } else{
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) SAA(s,a,i,j) = S_ya(i,j);
        S_ya = S_ya * get_S(P_ya, n_regions); //accumulate for next age
      }
    }
  }
  return SAA;
}
