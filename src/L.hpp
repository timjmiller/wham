template <class Type>
matrix<Type> get_L(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re) {
  /* 
     get "extra" mortality rate (n_years x n_regions)
       L_model: 0 = not used, 1 = constant (no RE), 2 = iid RE, 3 = AR1 RE (last two are distinquished on R side)
        L_pars: fixed effects parameters for extra mortality rate
            re: matrix of random effects to possibly use
  */
  int nr = L_model.size();
  int ny = re.rows();

  matrix<Type> L(ny,nr);
  L.setZero();
  for(int r = 0; r < n_regions; r++) if(L_model(r)>0) {
    for(int y = 0; y < ny; y++) {
      L(y,r) = exp(L_pars(r,0)); //mean
      if(L_model(r)>1){ //use random effects
        L(y,r) *= exp(re(y,r));
      }
    }
  }
  return(L);
}    
//done 

template <class Type>
vector<Type> get_nll_L(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re) {
  /* 
     get likelihood contributions for any RE used in "extra" mortality rate
       L_model: 0 = not used, 1 = constant (no RE), 2 = iid RE, 3 = AR1 RE (last two are distinquished on R side)
        L_pars: fixed effects parameters for extra mortality rate
            re: matrix of random effects to possibly use
  */
  int nr = L_model.size();
  int ny = re.rows();

  vector<Type> nll(nr);
  nll.setZero();
  for(int r = 0; r < n_regions; r++) if(L_model(r)>1) 
  {
    vector<Type> re_r = re.col(r);
    Type sig = exp(L_pars(r,1));
    Type rho = invlogit(L_pars(r,2),-1,1,1); //rho_trans(L_pars(r,2));
    nll(r) += SCALE(AR1(rho), sig)(re_r);
  }
  return(nll);
}
//done 

template <class Type>
matrix<Type> simulate_L_re(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re) {
  /* 
     simulate any RE used in "extra" mortality rate
       L_model: 0 = not used, 1 = constant (no RE), 2 = iid RE, 3 = AR1 RE (last two are distinquished on R side)
        L_pars: fixed effects parameters for extra mortality rate
            re: matrix of random effects to possibly use
  */
  int nr = L_model.size();
  int ny = re.rows();

  matrix<Type> re_sim = re;
  
  for(int r = 0; r < n_regions; r++) if(L_model(r)>1) 
  {
    vector<Type> re_r = re.col(r);
    Type sig = exp(L_pars(r,1));
    Type rho = invlogit(L_pars(r,2),-1,1,1); //rho_trans(L_pars(r,2));
    AR1(rho).simulate(re_r);
    re_r *= sig;
    for(int i = 0; i < re_r.size(); i++) re_sim(i,r) = re_r(i);
  }
  return(re_sim);
}
//done 
