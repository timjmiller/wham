template <class Type>
additive_ln_transform(vector<Type> x, int region, vector<int> do_move){
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
    } else { //prob of staying will be 1- prob of moving
      y(i) = 1.0;
    }
  }
  y /= sum(y);
  return(y);
}

template <class Type>
matrix<Type> get_P(int age, int year, int stock, int season, vector<int> fleet_regions, 
  matrix<int> can_move, int mig_type, Type time, array<Type> FAA, array<Type>MAA, 
  array<Type> trans_mu, matrix<Type> L)
//  vector<Type> mu, matrix<Type> L)
{
  /*Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2); */
  n_regions = MAA.dim(1);
  n_fleets = fleet_regions.size();
  
  int s = stock;
  int dim = n_regions+n_fleets+1;
  matrix<Type> P(dim,dim);
  P.setZero(); //zero it out.
  vector<Type> F(n_fleets), M(n_regions), Z(n_regions);
  
  for(int r = 0; r < n_regions; r++) 
  {
    M(r) = MAA(s,r,year,age) + L(year,r); //add in "extra" mortality
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
      matrix<Type> move(n_regions,n_regions);
      move.setZero();
      //additive logistic-normal transform for probability of movement out of regions only; prob staying is 1 - sum(prob move)
      //for(int i = 0; i < n_mu(stock,year,season,age); i++) move(mu_row(cum_n_mu + i)-1,mu_col(cum_n_mu + i)-1) = one/(one + exp(-mu(mu_pointer(cum_n_mu + i)-1)));
      
      for(int i = 0; i < n_regions; i++) {
        vector<Type> trans_par = trans_mu.row(i);
        vector<int> do_move = can_move.row(i);
        vector<Type> pmove = additive_transform(trans_par, i, do_move);
        for(int j = 0; j < n_regions; j++) move(i,j) = pmove(j);
      }

      if(time < 1e-15) 
      {
       //prob of survival is 1 when interval is zero
       for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) P(i,j) =  move(i,j);
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
            P(i,j) = exp(-Z(i) * time) * move(i,j); 
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
        for(int i = 0; i < n_regions; i++) A(i,dim-1) = M(i); //other dead
        for(int f = 0; f < n_fleets; f++) A(fleet_regions(f)-1,n_regions + f) = F(f);
        //for(int i = 0; i < n_mu(stock,year,season,age); i++) A(mu_row(cum_n_mu + i)-1,mu_col(cum_n_mu + i)-1) = exp(mu(mu_pointer(cum_n_mu + i)-1)); //log of transition intensities
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) if(can_move(i,j)==1){
          A(i,j) = exp(trans_mu(i,j)); //log of transition intensities
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
    } else
    {
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
  return P;
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

template <class Type>
matrix<Type> get_nll_N1(array<Type>log_N1, array<Type> N1_repars)
{
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  Type mu = N1_repars(s,r,0);
  Type sigma = exp(N1_repars(s,r,1));
  Type rho = rho_trans(N1_repars(s,r,2));
  nll_Ecov(i) += SCALE(AR1(Ecov_phi), Ecov_sig)(re_i);
  
}

template <class Type>
array<Type> set_N1(int N1_model, int n_ages, array<Type>Pmats)

  array<Type> N1(n_stocks, n_regions, n_ages);
  N1.setZero();
  if(N1_model == 0)
  {
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) 
    {
      NAA(s,r,0,a) = exp(log_N1_pars(s,r,a));
    }
  }
  /*if(N1_model == 1) //equilibrium approach
  {
    for(int r = 0; r < n_regions; r++) 
    {
      FAA_tot(r,0,0,a)
      for(int s = 0; s < n_stocks; s++) {
        NAA(s,r,0,0) = exp(log_N1_pars(0));
        for(int a = 1; a < n_ages; a++) 
        {
          matrix<Type> P = get_P(Pmats,s,r,0,a);
          if(a == n_ages-1) NAA(0,a) = NAA(s,r,0,a-1)/(1.0 + exp(-MAA(s,r,0,a) - exp(log_N1_pars(s,r,1)) * FAA_tot(0,a)/FAA_tot(0,n_ages-1)));
          else NAA(s,r,0,a) = NAA(s,r,0,a-1)* exp(-MAA(0,a) -  exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,n_ages-1));
        }
      }
    }
  } */
  if(N1_model == 2) //ar1 for each stock and region
