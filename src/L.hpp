template <class Type>
matrix<Type> get_L(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re, int n_years_model, int n_years_proj, 
    int proj_L_opt, vector<int> avg_years_ind) {
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
  for(int r = 0; r < nr; r++) if(L_model(r)>0) {
    for(int y = 0; y < n_years_model; y++) {
      L(y,r) = exp(L_pars(r,0)); //mean
      if(L_model(r)>1){ //use random effects
        L(y,r) *= exp(re(y,r));
      }
    }
  }
  if(n_years_proj > 0){
    if(proj_L_opt == 2){ // use average L over avg.yrs 
      vector<Type> L_proj = get_avg_L(L, avg_years_ind, 0);
      for(int y = n_years_model; y < ny; y++) for(int r = 0; r < nr; r++){
        L(y,r) = L_proj(r);
      }
    } else{
      for(int y = n_years_model; y < ny; y++) for(int r = 0; r < nr; r++) {
        L(y,r) = exp(L_pars(r,0)); //mean
        if(L_model(r)>1){ //use random effects
          L(y,r) *= exp(re(y,r));
        }
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
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int nr = L_model.size();

  vector<Type> nll(nr);
  nll.setZero();
  for(int r = 0; r < nr; r++) if(L_model(r)>1) 
  {
    vector<Type> re_r = re.col(r);
    Type mu = exp(L_pars(r,0));
    Type sig = exp(L_pars(r,1));
    Type rho = geninvlogit(L_pars(r,2),Type(-1),Type(1),Type(1)); //rho_trans(L_pars(r,2));
    nll(r) += SCALE(AR1(rho), sig)(re_r-mu);
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
  using namespace density; // necessary to use AR1, SCALE, SEPARABLE
  int nr = L_model.size();

  matrix<Type> re_sim = re;
  
  for(int r = 0; r < nr; r++) if(L_model(r)>1) 
  {
    vector<Type> re_r = re.col(r);
    Type mu = exp(L_pars(r,0));
    Type sig = exp(L_pars(r,1));
    Type rho = geninvlogit(L_pars(r,2),Type(-1),Type(1),Type(1)); //rho_trans(L_pars(r,2));
    AR1(rho).simulate(re_r);
    re_r *= sig;
    for(int i = 0; i < re_r.size(); i++) re_sim(i,r) = mu + re_r(i);
  }
  return(re_sim);
}
//done 
