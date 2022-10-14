/*
template <class Type>
matrix<Type> get_P(int age, int year, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, Type time, array<Type> FAA, array<Type>log_M, 
  array<Type> mu, matrix<Type> L) {
    /*
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
 /*
  int n_regions = log_M.dim(1);
  int n_fleets = fleet_regions.size();
  int dim = n_regions+n_fleets+1;
  matrix<Type> P(dim,dim);
  P.setZero(); //zero it out.
  vector<Type> F(n_fleets), M(n_regions), Z(n_regions);
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
          matrix<Type> A(dim,dim);
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

template <class Type>
matrix<Type> get_P(int age, int stock, int season, vector<int> fleet_regions, matrix<int> fleet_seasons,
  array<int> can_move, vector<int> mig_type, Type time, matrix<Type> FAA, array<Type>log_M, 
  array<Type> mu, vector<Type> L) {
    */
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
 /*
  array<Type> FAA_d(FAA.rows(),1, FAA.cols());
  FAA_d.setZero();
  for(int f = 0; f < FAA.rows(); f++) for(int a = 0; a < FAA.cols(); a++) FAA_d(f,0,a) = FAA(f,a);
  array<Type> log_M_d(log_M.dim(0),log_M.dim(1),1,log_M.dim(2));
  log_M_d.setZero();
  for(int s = 0; s < log_M.dim(0); s++) for(int a = 0; a < log_M.dim(2); a++) for(int r = 0; r < log_M.dim(1); r++) {
    log_M_d(s,r,0,a) = log_M(s,r,a);
  }
  array<Type> mu_d(mu.dim(0), mu.dim(1), mu.dim(2), 1, mu.dim(3),mu.dim(4));
  mu_d.setZero();
  for(int s = 0; s < mu.dim(0); s++) for(int a = 0; a < mu.dim(1); a++) for(int t = 0; t < mu.dim(2); t++) {
    for(int r = 0; r < mu.dim(3); r++) for(int rr = 0; rr < mu.dim(4); rr++){
      mu_d(s,a,t,0,r,rr) = mu_d(s,a,t,r,rr);
    }
  }
  matrix<Type> L_d(1,L.size());
  L_d.setZero();
  for(int r = 0; r < mu.dim(3); r++) L_d(0,r) = L(r);
  matrix<Type> P = get_P(age, 0, stock, season, fleet_regions, fleet_seasons, can_move, mig_type, time, FAA_d, log_M_d, mu_d, L_d);
  return(P);
}

template <class Type>
matrix<Type> get_S(matrix<Type> P, int n_regions){
  */
  /*
    extract the submatrix from a PTM that contains the proportions surviving in each region
              P: the probablity transition matrix
      n_regions: the number of regions
  */
 /*
  matrix<Type> S(n_regions,n_regions);
  S.setZero();
  for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) S(i,j) = P(i,j);
  return(S);
}
*/
template <class Type>
matrix<Type> get_D(matrix<Type> P, int n_regions, int n_fleets){
  /*
    extract the submatrix from a PTM that contains the proportions captured by each fleet in each region
              P: the probablity transition matrix
      n_regions: the number of regions
       n_fleets: the number of fleets
  */
  matrix<Type> D(n_regions,n_fleets);
  D.setZero();
  for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_fleets; j++) D(i,j) = P(i,j+n_regions);
  return(D);
}

template <class Type>
array<Type> get_annual_Ps(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the annual probability transition matrix for a given stock, age, season, year
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
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
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < n_seasons; t++) {
      //get numbers at age a for stock s in each region at time of spawning
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
    for(int i = 0; i < P_dim; i++) for(int j = 0; j < P_dim; j++) annual_Ps(s,y,a,i,j) = P_y(i,j);
  }
  return annual_Ps;
}

template <class Type>
array<Type> get_annual_SAA_spawn(vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, vector<Type> fracyr_seasons,
  matrix<Type> fracyr_SSB, vector<int> spawn_seasons, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the annual survival probabilities up to time of spawning for a given stock, age, season, year
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
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
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) {
    matrix<Type> P_y = I_mat; //reset for each year, age, stock
    for(int t = 0; t < spawn_seasons(s)-1; t++) { // only need to go up to the season prior to when spawning occurs
      //get numbers at age a for stock s in each region at time of spawning
      //P(t,u): PTM over entire season interval
      matrix<Type> P_t = get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
      P_y = P_y * P_t;
    }
    //P(0,t) x P(t_s-t): PTM over interval from to time of spawning within the season
    matrix<Type> P_SSB = P_y * get_P(a, y, s, spawn_seasons(s)-1, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_SSB(y,s), FAA, log_M, mu, L);
    for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) annual_SAA_SSB(s,y,a,i,j) = P_SSB(i,j);
  }
  return annual_SAA_SSB;
}

template <class Type>
array<Type> get_seasonal_Ps_y(int y, vector<int> fleet_regions, matrix<int> fleet_seasons, array<int> can_move, vector<int> mig_type, 
  vector<Type> fracyr_seasons, array<Type> FAA, array<Type> log_M, array<Type> mu, matrix<Type> L){
  /*
    produce the probability transition matrices for each stock, season, age for year y
      fleet_regions: n_fleets; which region each fleet is operating
      fleet_seasons: n_fleets x n_seasons; 0/1 indicating whether fleet is operating in the season
           can_move: n_stocks x ages x n_seasons x n_regions x n_regions; 0/1 determining whether movement can occur from one region to another
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
    matrix<Type> P_t = get_P(a, y, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons(t), FAA, log_M, mu, L);
    for(int d = 0; d < P_dim; d++) for(int dd = 0; dd < P_dim; dd++) P_seasonal_y(s,t,a,d,dd) = P_t(d,dd);
  }
  return P_seasonal_y;
}