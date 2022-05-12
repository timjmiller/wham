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
array<Type> get_NAA(int N1_model, array<Type> log_N1, array<Type>log_NAA, array<int> NAA_where){
  int n_stocks = log_N1.dim(0);
  int n_regions = log_N1.dim(1);
  int n_ages = log_N1.dim(3);
  int n_y = log_NAA.dim(2);
  array<Type> NAA(n_stocks, n_regions, n_y, n_ages);
  NAA.setZero();
  if(N1_model == 0 || N1_model == 2) {
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)==1){
      NAA(s,r,0,a) = exp(log_N1(s,r,a));
    }
  }
  for(int y = 1; y < n_y; y++){
    for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) if(NAA_where(s,r,a)==1){
      NAA(s,r,y,a) = exp(log_NAA(s,r,y-1,a));
    }
  }
  return(NAA);
}

template <class Type>
array<Type> get_SSB(array<Type>NAA_ssb, array<Type> waa, vector<int> waa_pointer_ssb, array<Type> mature){
  int n_stocks = NAA_ssb.dim(0);
  int n_ages = NAA_ssb.dim(2);
  int n_y = NAA_ssb.dim(1);
  
  matrix<Type> SSB(n_y, n_stocks) = get_SSB(NAA_ssb,waa,waa_pointer_ssb,mature);
  SSB.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_y; y++)
  {
    for(int a = 0; a < n_ages; a++) SSB(y,s) += NAA_SSB(s,y,a) * waa(waa_pointer_ssb(s)-1,y,a) * mature(s,y,a);
  }
  return(SSB);
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

//get SR_a, SR_b
template <class Type>
matrix<Type> get_SR_log_a(vector<int> recruit_model, matrix<Type> mean_rec_pars, array<Type> Ecov_lm, array<Type>Ecov_how, array<Type> Ecov_where)
  
  matrix<Type> log_SR_a(n_stocks,n_y);
  // calculate stock-recruit parameters (steepness, R0, a, b)
  if(recruit_model > 2) //BH or Ricker SR
  {
    log_SR_a.fill(mean_rec_pars(0));
    //log_SR_b.fill(mean_rec_pars(1));
    vector<Type> SR_h_tf(SR_h.size()); //different transformations for BH and Ricker
    if(recruit_model == 3) //BH stock recruit
    {
      /*if(use_steepness == 1)
      {
        SR_h.fill(0.2 + 0.8/(1+exp(-mean_rec_pars(0)))); //SR_a * SPR0/(4.0 + SR_a*SPR0);
        SR_R0.fill(exp(mean_rec_pars(1))); //(SR_a - 1/SPR0) / SR_b;
        log_SR_a = log(4 * SR_h/(exp(log_SPR0)*(1 - SR_h)));
        log_SR_b = log((5*SR_h - 1)/((1-SR_h)*SR_R0*exp(log_SPR0)));
      }
      else
      {*/
      //}
      for(int i=0; i < n_Ecov; i++){
        if(Ecov_where(i,0) == 1){ // if ecov i affects recruitment
          for(int y = 0; y < n_y; y++)
          {
            // (1) "controlling" = dens-indep mortality or (4) "masking" = metabolic/growth (decreases dR/dS)
            if((Ecov_how(i) == 1) | (Ecov_how(i) == 4))
            {
              log_SR_a(y) += Ecov_lm(i,0,y,0);
            }
            // (2) "limiting" = carrying capacity or (4) "masking" = metabolic/growth (decreases dR/dS)
            /*if((Ecov_how(i) == 2) | (Ecov_how(i) == 4))
            {
              log_SR_b(y) += Ecov_lm(i,0,y,0);
            }*/
          }
        }
      }
      /*if(use_steepness != 1)
      {
        SR_h = exp(log_SR_a) * exp(log_SPR0)/(4.0 + exp(log_SR_a + log_SPR0));
        SR_R0 = (exp(log_SR_a) - 1/exp(log_SPR0)) / exp(log_SR_b);
      }*/
      //SR_h_tf = log(SR_h - 0.2) - log(1 - SR_h);
    }
    if(recruit_model>3) //Ricker stock recruit
    {
      /*if(use_steepness == 1)
      {
        SR_h.fill(0.2 + exp(mean_rec_pars(0)));
        SR_R0.fill(exp(mean_rec_pars(1)));
        log_SR_a = 1.25*log(5*SR_h) - log_SPR0;
        log_SR_b = log(1.25*log(5*SR_h)/(SR_R0*exp(log_SPR0)));
      }*/
      for(int i=0; i < n_Ecov; i++){
        if(Ecov_where(i,0) == 1){ // if ecov i affects recruitment
          for(int y = 0; y < n_years_model + n_years_proj; y++)
          {
            if(Ecov_how(i) == 1) // "controlling" = dens-indep mortality
            {
              log_SR_a(y) += Ecov_lm(i,0,y,0);
            }
            if(Ecov_how(i) == 4) // "masking" = metabolic/growth (decreases dR/dS)
            { //NB: this is not identical to Iles and Beverton (1998), but their definition can give negative values of "b"
              log_SR_b(y) += 1.0 + Ecov_lm(i,0,y,0);
            }
          }
        }
      }
      if(use_steepness != 1)
      {
        SR_h = 0.2 * exp(0.8*log(exp(log_SR_a) * exp(log_SPR0)));
        SR_R0 = log(exp(log_SR_a + log_SPR0))/(exp(log_SR_b + log_SPR0));
      }
      SR_h_tf = log(SR_h - 0.2);
    }
    vector<Type> log_SR_R0 = log(SR_R0);
    if(do_post_samp.sum()==0){
      ADREPORT(log_SR_a);
      ADREPORT(log_SR_b);
      ADREPORT(SR_h_tf);
      ADREPORT(log_SR_R0);
    }
    REPORT(log_SR_a);
    REPORT(log_SR_b);
    REPORT(SR_h_tf);
    REPORT(log_SR_R0);
  }


////fix me: below
template <class Type>
vector<Type> get_pred_recruit_y(int y, vector<int> recruit_model, matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, matrix<Type> log_SR_a, 
  matrix<Type> log_SR_b, matrix<int> use_Ecov, matrix<int> Ecov_how, array<Type> Ecov_lm, vector<int> spawn_regions){

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
  int n_stocks = NAA.dim(0);
  vector<Type> pred_recruit(n_stocks);
  pred_recruit.setZero();
  for(int s = 0; s < n_stocks; s++){
    if(recruit_model(s) == 1) // random walk
    {
      pred_recruit(s) = NAA(s,r,y-1,0);
    } else
    {
      if(recruit_model(s) == 2) // random about mean
      {
        pred_recruit(s) = exp(mean_rec_pars(s,0));
        int nE = use_Ecov.dim(0); //rows
        for(int i=0; i < nE; i++){
          if(use_Ecov(i,s) == 1) if(Ecov_how(i,s) == 1) pred_recruit(s) *= exp(Ecov_lm(s,spawn_regions(s),y,i));
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
  return(pred_recruit);
}

template <class Type>
array<Type> get_pred_NAA(vector<int> recruit_model, matrix<Type> mean_rec_pars, matrix<Type> SSB, array<Type> NAA, matrix<Type> log_SR_a, 
  matrix<Type> log_SR_b, matrix<int> use_Ecov, matrix<int> Ecov_how, array<Type> Ecov_lm, vector<int> spawn_regions, array<Type> Ps){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: matrix of yearly SSB by stock (uses y-1 for any S-R relationship)
   * NAA: array of numbers at age
   * log_SR_a: matrix of yearly "a" parameters for SR function by stock 
   * log_SR_b: matrix of yearly "b" parameters for SR function by stock 
   * use_Ecov: matrix of 0/1 integer with first column determining if Ecov i affects recruitment (= 1)
   * Ecov_how: integer matrix with an element that tells how each Ecov is affecting recruitment for each stock
   * Ecov_lm: array that holds linear predictor for Ecov on recruitment
   * Ps: array of annual PTMs by stock and age
   */
  int n_stocks = NAA.dim(0);
  int n_years = NAA.dim(2);
  int n_ages = NAA.dim(3);
  int n_regions = NAA.dim(1);
  array<Type> pred_NAA(n_stocks,n_regions,n_years,n_ages);
  pred_NAA.setZero();

  // Expected recruitment
  for(int y = 1; y < n_years < y++){
    vector<Type> tmp = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
      log_SR_b, Ecov_where, Ecov_how, Ecov_lm, spawn_regions);
    for(int s = 0; s < n_stocks; s++){
      pred_NAA(s,spawn_regions(s),y,0) = tmp(s);
    }
    // calculate pred_NAA for ages after recruitment
    for(int a = 1; a < n_ages; a++) for(int s = 0; s < n_stocks; s++) {
      //vector<Type> NAA_s_last(n_regions);
      for(int r = 0; r < n_regions; r++) for(int r = 0; r < n_regions; r++){
      //pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
        pred_NAA(s,r,y,a) += Ps(s,y-1,a-1,rr,r) * NAA(s,y-1,a-1,rr);
      }
    }
    //plus group
    //pred_NAA(n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
    for(int r = 0; r < n_regions; r++) for(int r = 0; r < n_regions; r++){
        pred_NAA(s,r,y,n_ages-1) += Ps(s,y-1,n_ages-1,rr,r) * NAA(s,y-1,n_ages-1,rr);
    }
  }
  return(pred_NAA);
}


