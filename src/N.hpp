template <class Type>
matrix<Type> get_nll_N1(array<Type>log_N1, array<Type> N1_repars, matrix<int> N1_where)
{
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  matrix<Type> nll(s,r);
  nll.setZero();
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) if(N1_where(s,r)==1){
    Type mu = N1_repars(s,r,0);
    Type sigma = exp(N1_repars(s,r,1)); //marginal variance
    Type rho = rho_trans(N1_repars(s,r,2));
    vector<Type> re_sr(n_ages);
    for(int a = 0; a < n_ages; a++) re_sr(a) = log_N1(s,r,a) - mu
    nll(s,r) += SCALE(AR1(rho), sigma)(re_sr);
  }
}

array<Type> simulate_log_N1(array<Type>log_N1, array<Type> N1_repars, matrix<int> N1_where)
{
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  array<Type> sim_log_N1(s,r,a);
  sim_log_N1.setZero();
  for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) if(N1_where(s,r)==1){
    Type mu = N1_repars(s,r,0);
    Type sigma = exp(N1_repars(s,r,1)); //marginal variance
    Type rho = rho_trans(N1_repars(s,r,2));
    vector<Type> re_sr(n_ages);
    AR1(rho).simulate(re_sr);
    re_sr *= sigma;
    for(int a = 0; a < n_ages; a++) sim_log_N1(s,r,a) = re_sr(a) + mu
  }
  return(simulate_log_N1);
}

template <class Type>
array<Type> get_N1(int N1_model, array<Type> log_N1, matrix<int> N1_where){
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  array<Type> N1(n_stocks, n_regions, n_ages);
  N1.setZero();
  if(N1_model == 0 || N1_model == 2)
  {
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(N1_where(s,r)==1){
      N1(s,r,a) = exp(log_N1(s,r,a));
    }
  }
  return(N1);
}

template <class Type>
array<Type> get_pred_N1(int N1_model, array<Type> N1, matrix<int> N1_where){
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  array<Type> pred_N1(n_stocks, n_regions, n_ages);
  pred_N1.setZero();
  if(N1_model == 0)
  {
    pred_N1 = N1;
  } else {
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(N1_where(s,r)==1){
      pred_N1(s,r,a) = exp(N1_repars(s,r,0)); //exp of mean of AR1 process
    }
  }
  return(N1);
}

template <class Type>
matrix<Type> get_P(int age, int year, int stock, int season, vector<int> fleet_regions, 
  matrix<int> can_move, int mig_type, Type time, array<Type> FAA, array<Type>log_M_base, 
  array<Type> mu, matrix<Type> L)
{
  n_regions = log_M_base.dim(1);
  n_fleets = fleet_regions.size();
  int dim = n_regions+n_fleets+1;
  matrix<Type> P(dim,dim);
  P.setZero(); //zero it out.
  vector<Type> F(n_fleets), M(n_regions), Z(n_regions);
  
  for(int r = 0; r < n_regions; r++) 
  {
    M(r) = exp(log_M_base(stock,r,year,age)) + L(year,r); //add in "extra" mortality
    Z(r) = M(r);
  }
  for(int f = 0; f < n_fleets; f++) 
  {
    F(f) = FAA(year,season,f,age);
    Z(fleet_regions(f)-1) += F(f);
  }

  //if(n_mu(stock,year,season,age)>0) //migration is happening
  if(can_move.sum()>0) //migration is happening
  {
    if(mig_type == 0) //migration is instantaneous after survival and mortality, so P is easy.
    {
      if(time < 1e-15) 
      {
       //prob of survival is 1 when interval is zero
       for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) P(i,j) =  mu(stock,age,season,year,i,j);
      }
      else
      {
        for(int f = 0; f < n_fleets; f++) P(fleet_regions(f)-1,n_regions + f) = F(f) * (1.0 - exp(-Z(fleet_regions(f)-1) * time))/Z(fleet_regions(f)-1);
        for(int i = 0; i < n_regions; i++)
        {
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
    if(mig_type == 1) //migration occurs continuously during interval, so P is not so easy.
    {
       //prob of survival is 1 when interval is zero
      if(time < 1e-15) for(int i = 0; i < n_regions; i++) P(i,i) = 1.0;
      else
      {
        matrix<Type> A(dim,dim);
        A.setZero(); //zero it out.
        for(int i = 0; i < n_regions; i++) A(i,dim-1) += M(i); //other dead
        for(int f = 0; f < n_fleets; f++) A(fleet_regions(f)-1,n_regions + f) += F(f);
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) {
          //need to not put in mu diagonal because it has -sum of mig rates in there already.
          if(i != j) A(i,j) += mu(stock,age,season,year,i,j)); // transition intensities
        }
        for(int i = 0; i < n_regions; i++) A(i,i) = -(A.row(i)).sum(); //hazard
        A = A * time;
        P = expm(A);
      }
    }
  }
  else //no migration during this interval, so P is easy.
  {
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
  return(P);
}

//extract proportions surviving in each region over the interval.
template <class Type>
matrix<Type> get_S(matrix<Type> P, int n_regions){
  matrix<Type> S(n_regions,n_regions);
  for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++){
    S(i,j) = P(i,j);
  }
  return(S);
}

////fix me: below
template <class Type>
Type get_pred_recruit_y(int y, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSB, matrix<Type> NAA, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: matrix of 0/1 with first column determining if Ecov is affecting recruitment
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: array that holds linear predictor for Ecov
   */
  //recruit_model == 1, random walk
  Type pred_recruit = NAA(y-1,0);
  if(recruit_model == 1) // random walk
  {
    //pred_NAA(y,0) = NAA(y-1,0);
  }
  else
  {
    if(recruit_model == 2) // random about mean
    {
      pred_recruit = exp(mean_rec_pars(0));
      int nE = Ecov_where.rows();
      for(int i=0; i < nE; i++){
        if(Ecov_where(i,0) == 1) if(Ecov_how(i) == 1) pred_recruit *= exp(Ecov_lm(i,0,y,0));
      }
      //pred_NAA(y,0) = exp(mean_rec_pars(0));
      //if(Ecov_recruit > 0) if(Ecov_how(Ecov_recruit-1) == 1) pred_NAA(y,0) *= exp(Ecov_lm(y,Ecov_recruit-1));
    }
    else
    {
      if(recruit_model == 3) // BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
      {
        pred_recruit = exp(log_SR_a(y)) * SSB(y-1)/(1 + exp(log_SR_b(y))*SSB(y-1));
        //pred_NAA(y,0) = exp(log_SR_a(y)) * SSB(y-1)/(1 + exp(log_SR_b(y))*SSB(y-1));
      }
      else // recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
      {
        pred_recruit = exp(log_SR_a(y)) * SSB(y-1) * exp(-exp(log_SR_b(y)) * SSB(y-1));
        //pred_NAA(y,0) = exp(log_SR_a(y)) * SSB(y-1) * exp(-exp(log_SR_b(y)) * SSB(y-1));
      }
    }
  }
  return(pred_recruit);
}


