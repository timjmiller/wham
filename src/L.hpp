template <class Type>
matrix<Type> get_L(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re) {
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

template <class Type>
vector<Type> get_nll_L(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re)
{
  int nr = L_model.size();
  int ny = re.rows();

  vector<Type> nll(nr);
  nll.setZero();
  for(int r = 0; r < n_regions; r++) if(L_model(r)>1) 
  {
    vector<Type> re_r = re.col(r);
    Type sig = exp(L_pars(r,1));
    Type rho = rho_trans(L_pars(r,2));
    nll(r) += SCALE(AR1(rho), sig)(re_r);
  }
  return(nll);
}

template <class Type>
matrix<Type> simulate_L_re(vector<int> L_model, matrix<Type> L_pars, matrix<Type> re) {
  int nr = L_model.size();
  int ny = re.rows();

  matrix<Type> re_sim = re;
  
  for(int r = 0; r < n_regions; r++) if(L_model(r)>1) 
  {
    vector<Type> re_r = re.col(r);
    Type sig = exp(L_pars(r,1));
    Type rho = rho_trans(L_pars(r,2));
    AR1(rho).simulate(re_r);
    re_r *= sig;
    for(int i = 0; i < re_r.size(); i++) re_sim(i,r) = re_r(i);
  }
  return(re_sim);
}
