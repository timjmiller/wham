template <class Type>
vector<Type> get_nll_sel(vector<int> selblock_models_re, vector<int> n_years_selblocks, vector<int> n_selpars_est, 
  array<Type> selpars_re, matrix<Type> sel_repars){
  /* 
     get nll contribtions for any selectivity random effects
        selblock_models_re: (n_selblocks) 1 = no RE, 2 = IID, 3 = ar1, 4 = ar1_y, 5 = 2dar1
         n_years_selblocks: for each block, number of years the block covers
             n_selpars_est: n_selbocks, how many selpars are actually estimated (not fixed at 0 or 1) 
                selpars_re: (n_selbocks x n_years x n_ages) deviations in selectivity parameters (random effects), length = sum(n_selpars)*n_years per block
                sel_repars: parameters controlling selpars_re, dim = n_blocks, 3 (sigma, rho, rho_y)
  */
  using namespace density;
  int n_selblocks = selblock_models_re.size();
  vector<Type> nll_sel(n_selblocks);
  nll_sel.setZero();
  //int istart = 0;
  for(int b = 0; b < n_selblocks; b++){

    if(selblock_models_re(b) > 1){
      // fill in sel devs from RE vector, selpars_re (fixed at 0 if RE off)
      array<Type> tmp(n_years_selblocks(b), n_selpars_est(b));
      for(int i = 0; i < n_years_selblocks(b); i++) for(int j=0; j<n_selpars_est(b); j++){
        tmp(i,j) = selpars_re(b,i,j);
        //tmp.col(j) = selpars_re.segment(istart,n_years_selblocks(b));
        //istart += n_years_selblocks(b);
      }

      //question: is it faster here to just work on the selectivity parameters as re rather than the deviations?
      // likelihood of RE sel devs (if turned on)
      Type sigma = exp(sel_repars(b,0)); // sd selectivity deviations (fixed effect)
      //rho_trans ensures correlation parameter is between -1 and 1, see helper_functions.hpp
      Type rho = geninvlogit(sel_repars(b,1),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy // among-par correlation selectivity deviations (fixed effect) 
      Type rho_y = geninvlogit(sel_repars(b,2),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy // among-year correlation selectivity deviations (fixed effect)
      
      if((selblock_models_re(b) == 2) | (selblock_models_re(b) == 5)){
        // 2D AR1 process on selectivity parameter deviations
        Type Sigma_sig_sel = pow(pow(sigma,2) / ((1-pow(rho_y,2))*(1-pow(rho,2))),0.5);
        nll_sel(b) += SCALE(SEPARABLE(AR1(rho),AR1(rho_y)), Sigma_sig_sel)(tmp);
      } else {
        // 1D AR1 process on selectivity parameter deviations
        if(selblock_models_re(b) == 3){ // ar1 across parameters in selblock, useful for age-specific pars.
          vector<Type> tmp0 = tmp.matrix().row(0); //random effects are constant across years 
          Type Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho,2)),0.5);
          nll_sel(b) += SCALE(AR1(rho), Sigma_sig_sel)(tmp0);
        } else { // selblock_models_re(b) = 4, ar1_y, not sure if this one really makes sense.
          vector<Type> tmp0 = tmp.matrix().col(0); //random effects are constant within years 
          Type Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho_y,2)),0.5);
          //Sigma_sig_sel = sigma;
          nll_sel(b) += SCALE(AR1(rho_y), Sigma_sig_sel)(tmp0);
        }
      }
    }
  }
  return nll_sel;
}

template <class Type>
array<Type> simulate_selpars_re(vector<int> selblock_models_re, vector<int> n_years_selblocks, vector<int> n_selpars_est, 
  array<Type> selpars_re, matrix<Type> sel_repars){
  /* 
     simulate any selectivity random effects
        selblock_models_re: (n_selblocks) 1 = no RE, 2 = IID, 3 = ar1, 4 = ar1_y, 5 = 2dar1
         n_years_selblocks: for each block, number of years the block covers
             n_selpars_est: n_selbocks, how many selpars are actually estimated (not fixed at 0 or 1) 
                selpars_re: (n_selbocks x n_years x n_ages) deviations in selectivity parameters (random effects), length = sum(n_selpars)*n_years per block
                sel_repars: parameters controlling selpars_re, dim = n_blocks, 3 (sigma, rho, rho_y)
  */
  using namespace density;
  int n_selblocks = selblock_models_re.size();
  //int istart = 0;
  array<Type> sim_selpars_re = selpars_re;
  for(int b = 0; b < n_selblocks; b++){

    if(selblock_models_re(b) > 1){
      // fill in sel devs from RE vector, selpars_re (fixed at 0 if RE off)
      array<Type> tmp(n_years_selblocks(b), n_selpars_est(b));
      for(int i = 0; i < n_years_selblocks(b); i++) for(int j=0; j<n_selpars_est(b); j++){
        tmp(i,j) = selpars_re(b,i,j);
        //tmp.col(j) = selpars_re.segment(istart,n_years_selblocks(b));
        //istart += n_years_selblocks(b);
      }

      //question: is it faster here to just work on the selectivity parameters as re rather than the deviations?
      // likelihood of RE sel devs (if turned on)
      Type sigma = exp(sel_repars(b,0)); // sd selectivity deviations (fixed effect)
      //rho_trans ensures correlation parameter is between -1 and 1, see helper_functions.hpp
      Type rho = geninvlogit(sel_repars(b,1),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy // among-par correlation selectivity deviations (fixed effect) 
      Type rho_y = geninvlogit(sel_repars(b,2),Type(-1),Type(1),Type(1));//using scale =1 ,2 is legacy // among-year correlation selectivity deviations (fixed effect)
      Type Sigma_sig_sel = 0;

      if((selblock_models_re(b) == 2) | (selblock_models_re(b) == 5)){
        // 2D AR1 process on selectivity parameter deviations
        Sigma_sig_sel = pow(pow(sigma,2) / ((1-pow(rho_y,2))*(1-pow(rho,2))),0.5);
        SEPARABLE(AR1(rho),AR1(rho_y)).simulate(tmp);
      } else {
        // 1D AR1 process on selectivity parameter deviations
        if(selblock_models_re(b) == 3){ // ar1 across parameters in selblock, useful for age-specific pars.
          vector<Type> tmp0 = tmp.matrix().row(0); //random effects are constant across years 
          Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho,2)),0.5);
          AR1(rho).simulate(tmp0);
          for(int y = 0; y < tmp.rows(); y++) for(int i = 0; i < tmp0.size(); i++) tmp(y,i) = tmp0(i);
        } else { // selblock_models_re(b) = 4, ar1_y, not sure if this one really makes sense.
          vector<Type> tmp0 = tmp.matrix().col(0); //random effects are constant within years 
          Sigma_sig_sel = pow(pow(sigma,2) / (1-pow(rho_y,2)),0.5);
          AR1(rho_y).simulate(tmp0);
          for(int a = 0; a < tmp.cols(); a++) tmp.col(a) = tmp0;
        }
      }
      tmp = tmp * Sigma_sig_sel;
      //istart -= n_selpars_est(b) * n_years_selblocks(b); //bring it back to the beginning for this selblock
      for(int j=0; j<n_selpars_est(b); j++){
        for(int y = 0; y < n_years_selblocks(b); y++){
          sim_selpars_re(b,y,j) = tmp(y,j);
          //sim_selpars_re(istart) = tmp(y,j);
          //istart++;
        }
      }
    }
  }
  return sim_selpars_re;
}

template <class Type>
vector<matrix<Type> > get_selpars_re_mats(vector<int> n_selpars, matrix<int> selblock_years, matrix<int> selpars_est, 
  int n_years_model, array<Type> selpars_re, vector<int> selblock_models, vector<int> selblock_models_re){
  /* 
    get vector of matrices of selectivity random effects.
                 n_selpars: n_selblocks. how many mean selectivity parameters estimated for each selblock 
            selblock_years: n_years_model x n_selblocks, = 1 if block covers year, = 0 if not
               selpars_est: n_blocks x (n_pars(6) + n_ages), 0/1; is the selpar estimated in this block?
             n_years_model: number of non-projection years in the model
                selpars_re: (n_selbocks x n_years x n_ages) deviations in selectivity parameters (random effects), length = sum(n_selpars)*n_years per block
           selblock_models: n_selblocks. which (mean) selectivity model for each block
        selblock_models_re: (n_selblocks) 1 = no RE, 2 = IID, 3 = ar1, 4 = ar1_y, 5 = 2dar1
  */
  
  int n_selblocks = n_selpars.size();
  int n_ages = selpars_est.cols() - 6;
  vector<matrix<Type> > selpars_re_mats(n_selblocks);
  //int istart = 0;
  int ct = 0;
  for(int b = 0; b < n_selblocks; b++){
    matrix<Type> tmp2(n_years_model, n_selpars(b));
    tmp2.setZero();
    selpars_re_mats(b) = tmp2;

    int jstart = 0; // offset for indexing selectivity pars, depends on selectivity model for block b: n_ages (age-specific) + 2 (logistic +/-) + 4 (double-logistic)
    if((selblock_models(b) == 2) | (selblock_models(b) == 4)) jstart = n_ages;
    if(selblock_models(b) == 3) jstart = n_ages + 2;

    if(selblock_models_re(b) > 1){
      // construct deviations array with full dimensions (n_years_model instead of n_years_selblocks, n_selpars instead of n_selpars_est)
      int jj = 0;
      for(int j=0; j<n_selpars(b); j++){
        if(selpars_est(b,j+jstart) > 0){
          ct = 0;
          for(int y=0; y<n_years_model; y++){
            if(selblock_years(y,b) == 1){
              selpars_re_mats(b)(y,j) = selpars_re(b,ct,jj);
              ct++;
            }
          }
          jj++;
        }
      }
    }
  }
  return selpars_re_mats; //even if not simulated
}

template <class Type>
vector<matrix<Type> > get_selpars(vector<int> selblock_models, vector<int> n_selpars, matrix<Type> logit_selpars, 
  vector<matrix<Type> >  selpars_re_mats, matrix<Type> selpars_lower, matrix<Type> selpars_upper, int n_years_model){
  /* 
    get vector of matrices of selectivity parameters.
      selblock_models: n_selblocks. which (mean) selectivity model for each block
            n_selpars: n_selblocks. how many selectivity parameters for each selblock 
        logit_selpars: n_selblocks x (6+n_ages) matrix of mean logit-selectivity parameters
      selpars_re_mats: vector of matrices of selectivity random effects
        selpars_lower: n_selblocks x (6+n_ages) lower bound of selectivity parameters for invlogit transformation (default = 0)
        selpars_upper: n_selblocks x (6+n_ages) upper bound of selectivity parameters for invlogit transformation (default = 1 or n_ages)
        n_years_model: number of non-projection years in the model
  */

  int n_selblocks = selblock_models.size();
  int n_ages = logit_selpars.cols() - 6;
  vector<matrix<Type> > selpars(n_selblocks); // selectivity parameter matrices for each block, nyears x npars
  for(int b = 0; b < n_selblocks; b++){
    int jstart = 0; // offset for indexing selectivity pars, depends on selectivity model for block b: n_ages (age-specific) + 2 (logistic) + 4 (double-logistic)
    if((selblock_models(b) == 2) | (selblock_models(b) == 4)) jstart = n_ages;
    if(selblock_models(b) == 3) jstart = n_ages + 2;

    // get selpars = mean + deviations
    matrix<Type> tmp1(n_years_model, n_selpars(b));
    for(int j=jstart; j<(jstart+n_selpars(b)); j++){ // transform from logit-scale
      for(int i=0; i<n_years_model; i++){
        Type logit_sel_re = logit_selpars(b,j) + selpars_re_mats(b)(i,j-jstart);
        tmp1(i,j-jstart) = geninvlogit(logit_sel_re,selpars_lower(b,j), selpars_upper(b,j),Type(1));
        //tmp1(i,j-jstart) = selpars_lower(b,j) + (selpars_upper(b,j) - selpars_lower(b,j)) / (1.0 + exp(-(logit_selpars(b,j) + selpars_re_mats(b).matrix()(i,j-jstart))));
      }
    }
    selpars(b) = tmp1;
  }
  return selpars;

}

template <class Type>
vector<matrix<Type> > get_selAA(int n_years, int n_ages, int n_selblocks, vector<matrix<Type> > selpars, 
  vector<int> selblock_models) {
  /* 
    get vector of matrices of selectivity at age.
              n_years: n_years_model 
               n_ages: n_ages
          n_selblocks: n_selblocks
              selpars: vector of matrices of selectivity parameters
      selblock_models: n_selblocks. which (mean) selectivity model for each block
  */
  vector<matrix<Type> > selAA(n_selblocks);
  for(int b = 0; b < n_selblocks; b++)
  {
    matrix<Type> tmp(n_years, n_ages);
    if(selblock_models(b) == 1) tmp = selpars(b); //proportions at age
    else
    { //logistic or double-logistic
      if(selblock_models(b) == 2)
      { //increasing logistic
        for(int y = 0; y < n_years; y++)
        {
          Type a50 = selpars(b)(y,0); // a50 parameter in year y
          Type k = selpars(b)(y,1); //  1/slope in year y
          Type age = 0.0;
          for(int a = 0; a < n_ages; a++)
          {
            age += 1.0;
            tmp(y,a) = 1.0/(1.0 + exp(-(age - a50)/k));
          }
          for(int a = 0; a < n_ages; a++) tmp(y,a) = tmp(y,a)/tmp(y,n_ages-1);
        }
      }
      else
      { //double logistic
        if(selblock_models(b) == 3)
        {
          for(int y = 0; y < n_years; y++)
          {
            Type a50_1 = selpars(b)(y,0); // a50 parameter in year y
            Type k_1 = selpars(b)(y,1); //  1/slope in year y
            Type a50_2 = selpars(b)(y,2);
            Type k_2 = selpars(b)(y,3);
            Type age = 0.0;
            for (int a = 0; a < n_ages; a++)
            {
              age += 1.0;
     	        tmp(y,a) = 1.0/(1.0 + exp(-(age - a50_1)/k_1));
              tmp(y,a) *= 1.0/(1.0 + exp((age - a50_2)/k_2)); //1-p
            }
          }
        }
        else //model 4: declining logistic
        {
          for(int y = 0; y < n_years; y++)
          {
            Type a50 = selpars(b)(y,0); // a50 parameter in year y
            Type k = selpars(b)(y,1); //  1/slope in year y
            Type age = 0.0;
            for (int a = 0; a < n_ages; a++)
            {
              age += 1.0;
              tmp(y,a) = 1.0/(1.0 + exp((age - a50)/k));
            }
            for (int a = 0; a < n_ages; a++) tmp(y,a) = tmp(y,a)/tmp(y,0);
          }
        }
      }
    }
    selAA(b) = tmp;
  }
  return selAA;
}
