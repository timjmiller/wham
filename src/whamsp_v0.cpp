#include <TMB.hpp>
#include <iostream>
#include "helper_functions.hpp"
#include "age_comp_osa.hpp"
#include "age_comp_sim.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_years_model); 
  DATA_INTEGER(n_seasons);
  DATA_VECTOR(fracyr_seasons); //length of intervals for seasons
  DATA_INTEGER(n_regions);
  DATA_INTEGER(n_stocks);
  DATA_INTEGER(n_ages_model);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_IVECTOR(mig_type); //n_stocks. 0 = migration after survival, 1 = movement and mortality simultaneous
  
  DATA_MATRIX(fracseason_SSB); //n_years x n_stocks
  DATA_MATRIX(fracyr_SSB); //n_years x n_stocks
  DATA_IVECTOR(spawn_regions); //n_stocks
  DATA_IVECTOR(spawn_seasons); //n_stocks
  DATA_IARRAY(n_mu) //n_stocks x n_years x x n_seasons x n_ages
  DATA_IVECTOR(mu_row) //length the total number of migration parameters
  DATA_IVECTOR(mu_col) //length the total number of migration parameters
  DATA_IVECTOR(mu_pointer) //n_stocks * n_years * n_ages
  DATA_INTEGER(n_seasons_recruited); //easiest if this is = to n_seasons 
  DATA_ARRAY(mature); //n_stocks x n_years x n_ages
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_IVECTOR(waa_pointer_totcatch); //n_regions
  DATA_IVECTOR(waa_pointer_indices);
  DATA_IVECTOR(waa_pointer_ssb); //n_stocks
  DATA_ARRAY(waa); //(n_fleets + n_indices + n_stocks + 1(totcatch)) x n_years x n_ages_model
  //M is a parameter
  //DATA_ARRAY(MAA); // n_years x n_ages x n_regions
  
  DATA_IVECTOR(fleet_regions); //length = n_fleets
  DATA_IMATRIX(fleet_seasons); //n_fleets x n_seasons; indicator with length n_seasons, e.g., 1,1,1.. would mean operating all year.
  DATA_MATRIX(agg_catch); //n_years x n_fleets 
  DATA_IMATRIX(use_agg_catch); //n_years x n_fleets 
  DATA_MATRIX(agg_catch_sigma); //n_years x n_fleets 
  DATA_ARRAY(catch_paa); //n_fleets x n_years x n_ages_model 
  DATA_IMATRIX(use_catch_paa); //n_years x n_fleets 
  DATA_MATRIX(catch_Neff); //n_years x n_fleets 
  DATA_IVECTOR(age_comp_model_fleets); //length = n_fleets
  
  DATA_IVECTOR(index_regions); //length = n_indices
  DATA_IVECTOR(index_seasons); //n_indices x n_seasons. indicator must only be for one season
  DATA_IVECTOR(units_indices); //length = n_indices
  DATA_MATRIX(fracyr_indices); //n_years x n_indices
  DATA_MATRIX(agg_indices); //n_years x n_indices
  DATA_IMATRIX(use_indices);  //n_years x n_indices
  DATA_MATRIX(agg_index_sigma); //n_years x n_indices
  DATA_IVECTOR(units_index_paa); //length = n_indices
  DATA_ARRAY(index_paa); //n_indices x n_years x n_ages_model
  DATA_IMATRIX(use_index_paa); //n_years x n_indices
  DATA_MATRIX(index_Neff); //n_years x n_indices
  DATA_IVECTOR(age_comp_model_indices); //length = n_indices
  DATA_VECTOR(q_lower); //length = n_indices
  DATA_VECTOR(q_upper); //length = n_indices
  DATA_IVECTOR(use_q_prior); //length = n_indices
  DATA_VECTOR(logit_q_prior_sigma); //length = n_indices
  DATA_IVECTOR(use_q_re);  //length = n_indices, 0= no re, >0 = use re   
  
  DATA_IVECTOR(selblock_models);
  int n_selblocks = selblock_models.size();
  DATA_IVECTOR(selblock_models_re); // for each block: 1 = none, 2 = IID, 3 = ar1, 4 = ar1_y, 5 = 2dar1
  DATA_IVECTOR(n_selpars); //length = n_selblocks
  DATA_IMATRIX(selpars_est); // n_blocks x (n_pars + n_ages), is the selpar estimated in this block?
  DATA_IVECTOR(n_selpars_est); // of the selpars, how many are actually estimated (not fixed at 0 or 1)
  DATA_IVECTOR(n_years_selblocks); // for each block, number of years the block covers
  DATA_IMATRIX(selblock_years); // n_years_model x n_selblocks, = 1 if block covers year, = 0 if not
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_MATRIX(selpars_lower);
  DATA_MATRIX(selpars_upper);
  
  DATA_IVECTOR(recruit_model); //length = n_stocks
  DATA_IVECTOR(n_M_a); //length = n_stocks
  DATA_IVECTOR(M_model); //length = n_stocks; 1: "constant", 2: "age-specific", 3: "weight-at-age"
  DATA_IVECTOR(N1_model); //length = n_stocks; 0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations
  DATA_IVECTOR(M_re_model); //length = n_stocks; 1 = none, 2 = IID, 3 = ar1_a, 4 = ar1_y, 5 = 2dar1
  DATA_IVECTOR(use_b_prior); //length = n_stocks; for M_model = 3: M = a*W^b model

  DATA_IMATRIX(which_F_age); // (n_years_model + n_years_proj) x n_stocks; which age of F to use for full total F for msy/ypr calculations and projections
  DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  
  //DATA_IMATRIX(R1_pointer); //n_stocks x n_regions. Tells where log_R1 parameters are used.
  //DATA_IMATRIX(N1_sigma_pointer); //n_stocks x n_regions. Tells which N1_sigma_par to use where.
  //DATA_IARRAY(NAA_sigma_pointers); //n_stocks x n_ages x n_regions
  DATA_IARRAY(NAA_re_indicator); //n_stocks x (n_years-1) x n_ages x n_regions will estimate a random effect where indicator is not 0. sum = n_NAA_re.
  
  PARAMETER_MATRIX(mean_rec_pars); //n_stocks x n_rec_pars (determined by recruit_model)
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(log_F1);
  PARAMETER_MATRIX(F_devs);
  PARAMETER_VECTOR(mu); //migration parameters
  //N1 might need some tweeking. for example, if all fish are forced to be in spawning region at the beginning of the year, then there should be no N1 in other regions.
  PARAMETER_VECTOR(log_R1); //length must be consistent with R1_pointer above.
  PARAMETER_VECTOR(log_N1_sigma); //length must be consistent with N1_sigma_pointer above. Used if N1_model = 2 and log_N1 is a random walk in age.
  PARAMETER_VECTOR(log_NAA_sigma); // vector sigmas used with NAA_sigma_pointers
  PARAMETER_VECTOR(estimated_selpars);
  PARAMETER_VECTOR(catch_paa_pars);
  PARAMETER_VECTOR(index_paa_pars);
  //Depending on what we do here; the steps through the seasons will only be using NAA from the beginning of the year or NAA will vary through the seasons too
  //Just have annual NAA right now
  //need to use vector because NAA may be zero in some regions due to migration parameterization, need to figure out length in R
  PARAMETER_VECTOR(log_N1_re); //at most, length is n_stocks * n_regions * (n_ages -1) 
  PARAMETER_VECTOR(log_NAA_re); //n_NAA_re. depends on whether NAA in some regions are zero. dimension is n_NAA_re. used with NAA_re_indicator to fill NAA.
  //PARAMETER_MATRIX(log_NAA); //n_seasons x n_stocks* n_ages
  //std::cout << "here 1"  << "\n";
  PARAMETER_VECTOR(log_L1); //n_regions
  PARAMETER_MATRIX(log_L_re); //n_years-1 x n_regions
  PARAMETER_VECTOR(log_L_sigma); //n_regions
  
  int P_dim = n_regions + n_fleets + 1;
  matrix<Type> I_mat(P_dim,P_dim);
  for(int i = 0; i < P_dim; i++) I_mat(i,i) = 1.0;

  //Type zero = Type(0);
  //Type one = Type(1);
  //Type half = Type(0.5);
  //Type two = Type(2);
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }
  vector<Type> sigma2_log_NAA = exp(log_NAA_sigma*two);
  matrix<Type> SSB(n_years,n_stocks);
  matrix<Type> log_SSB(n_years,n_stocks);
  matrix<Type> F(n_years,n_fleets);
  matrix<Type> log_F(n_years,n_fleets);
  array<Type> pred_CAA(n_years,n_fleets,n_ages);
  array<Type> pred_stock_CAA(n_years,n_stocks,n_ages,n_regions);
  array<Type> pred_stock_catch(n_years,n_stocks,n_regions);
  matrix<Type> pred_catch_region(n_years,n_regions);
  array<Type> pred_stock_prop_catch(n_years,n_stocks,n_regions);
  array<Type> pred_catch_paa(n_years,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years,n_fleets);
  array<Type> pred_IAA(n_years,n_indices,n_ages);
  array<Type> pred_index_paa(n_years,n_indices,n_ages);
  matrix<Type> pred_indices(n_years,n_indices);
  matrix<Type> log_pred_catch(n_years,n_fleets);
  array<Type> NAA(n_stocks,n_years,n_ages,n_regions);
  array<Type> pred_NAA(n_stocks,n_years,n_ages,n_regions);
  array<Type> NAA_SSB(n_stocks,n_years,n_ages);
  array<Type> FAA(n_years,n_seasons,n_fleets,n_ages);
  array<Type> FAA_tot(n_years,n_seasons,n_ages,n_regions);
  array<Type> all_P(n_stocks,n_years,n_ages,P_dim,P_dim);
  array<Type> all_P_SSB(n_stocks,n_years,n_ages,n_regions,n_regions);
  array<Type> P_season_terminal(n_seasons,n_stocks,n_ages,P_dim,P_dim);
  array<Type> NAA_index(n_stocks,n_indices,n_years,n_ages);
  array<Type> QAA(n_years,n_indices,n_ages);
  matrix<Type> selblocks(n_selblocks,n_ages);
  vector<Type> q(n_indices);
  matrix<Type> L(n_years,n_regions);
  Type nll = zero; //negative log-likelihood
  vector<Type> t_paa(n_ages), t_pred_paa(n_ages);

  selblocks = get_selblocks(n_ages, n_selblocks, n_estimated_selpars, n_other_selpars, selblock_models, estimated_selpar_pointers, other_selpar_pointers, estimated_selpars, other_selpars, selpars_lower, selpars_upper);
  for(int i = 0; i < n_indices; i++)
  {
    q(i) = q_lower(i) + (q_upper(i) - q_lower(i))/(1 + exp(-logit_q(i)));
    for(int y = 0; y < n_years; y++) 
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(i) * selblocks(selblock_pointer_indices(y,i)-1,a);
    }
  }
  //std::cout << "here 2"  << "\n";
 
  L.setZero();
  Type nll_L = zero;
  for(int y = 0; y < n_years; y++) for(int r = 0; r < n_regions; r++) 
  {
    if(y ==0) L(y,r) = exp(log_L1(r));
    else 
    {
      L(y,r) = exp(log_L_re(y-1,r));
      nll_L -= dnorm(log_L_re(y-1,r), log(L(y-1,r)), exp(log_L_sigma(r)), 1);
    }
  }
  
  FAA_tot.setZero();
  FAA.setZero();
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
    for(int a = 0; a < n_ages; a++) 
    {
      
      for(int t = 0; t < n_seasons; t++) if(fleet_seasons(f,t) == 1) 
      {
        FAA(0,t,f,a) = F(0,f) * selblocks(selblock_pointer_fleets(0,f)-1,a);
        FAA_tot(0,t,a,fleet_regions(f)-1) = FAA_tot(0,t,a,fleet_regions(f)-1) + FAA(0,t,f,a);
      }
    }
    for(int y = 1; y < n_years; y++) 
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
      for(int a = 0; a < n_ages; a++) 
      {
        for(int t = 0; t < n_seasons; t++) if(fleet_seasons(f,t) == 1) 
        {
          FAA(y,t,f,a) = F(y,f) * selblocks(selblock_pointer_fleets(y,f)-1,a);
          FAA_tot(y,t,a,fleet_regions(f)-1) = FAA_tot(y,t,a,fleet_regions(f)-1) + FAA(y,t,f,a);
        }
      }
    }
  }
  //std::cout << "here 3"  << "\n";
  NAA.setZero();
  pred_NAA.setZero();
  Type nll_N1 = zero;
  /*for(int s = 0; s < n_stocks; s++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) 
  {
    if(N1_pointer(s,a,r) > 0) NAA(s,0,a,r) = exp(log_N1(N1_pointer(s,a,r)-1));
    pred_NAA(s,0,a,r) = NAA(s,0,a,r);
  }*/
  int N1_re_ind = 0;
  for(int s = 0; s < n_stocks; s++)  for(int r = 0; r < n_regions; r++) if(R1_pointer(s,r) > 0)
  {
    NAA(s,0,0,r) = pred_NAA(s,0,0,r) = exp(log_R1(R1_pointer(s,r)-1));
    for(int a = 1; a < n_ages; a++)  
    {
      NAA(s,0,a,r) = pred_NAA(s,0,a,r) = exp(log_N1_re(N1_re_ind));
      nll_N1 -= dnorm(log(NAA(s,0,a,r)), log(NAA(s,0,a-1,r)), exp(log_N1_sigma(N1_sigma_pointer(s,r)-1)),1);
      N1_re_ind ++;
    }
  }
  //for(int a = 0; a < n_ages; a++) see(NAA(0,0,a,0));
  //for(int a = 0; a < n_ages; a++) see(NAA(1,0,a,1));
  see(N1_re_ind);
  see(nll_N1);
  nll += nll_N1;
  int ind = 0;
  for(int s = 0; s < n_stocks; s++) for(int y = 1; y < n_years; y++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) 
  {
    if(NAA_re_indicator(s,y-1,a,r))
    {
      NAA(s,y,a,r) = exp(log_NAA_re(ind)); //random effects NAA
      ind ++;
    }
  }
  
  //std::cout << "here 4"  << "\n";
  int cum_n_mu = 0;
  
  //get probability transition matrix, NAA at spawning and NAA for each index.
  matrix<Type> P1(P_dim,P_dim), P_SSB(P_dim,P_dim), P_index(P_dim,P_dim);
  NAA_SSB.setZero();
  NAA_index.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++) for(int a = 0; a < n_ages; a++) 
  {
    P1 = I_mat;
    for(int t = 0; t < n_seasons; t++) 
    {
      //get numbers at age a for stock s in each region at time of spawning
      if(t == spawn_seasons(s)-1)
      {
        P_SSB = P1 * get_P(n_regions, n_fleets, a, y, s, t, fleet_regions, cum_n_mu, n_mu, mu_row, mu_col, mu_pointer, mig_type(s), fracyr_SSB(y,s), FAA, MAA, mu, L);
        for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) all_P_SSB(s,y,a,i,j) = P_SSB(i,j);
        //spawners only spawning in stock region...
        for(int r = 0; r < n_regions; r++) NAA_SSB(s,y,a) += P_SSB(r,spawn_regions(s)-1) * NAA(s,y,a,r);
      }
      for(int i = 0; i < n_indices; i++) 
      {
        if(t == index_seasons(i)-1)
        { 
          P_index = P1 * get_P(n_regions, n_fleets, a, y, s, t, fleet_regions, cum_n_mu, n_mu, mu_row, mu_col, mu_pointer, mig_type(s), fracyr_indices(y,i), FAA, MAA, mu, L);
          for(int r = 0; r < n_regions; r++) NAA_index(s,i,y,a) += P_index(r,index_regions(i)-1) * NAA(s,y,a,r);
          /*if(y == 0 && i == 1 && a == 0)
          {
            std::cout << "P_index for index " << i << " in year " << y  << "\n";
            std::cout << P_index << "\n";
            std::cout << index_regions(i) << "\n";
            std::cout << NAA_index(s,i,y,a) << "\n";
            std::cout << NAA(s,y,a,0) << "\n";
          } */
        }
      }
      matrix<Type> P_t = get_P(n_regions, n_fleets, a, y, s, t, fleet_regions, cum_n_mu, n_mu, mu_row, mu_col, mu_pointer, mig_type(s), fracyr_seasons(t), FAA, MAA, mu, L);
      if(y == n_years-1) 
      {
        for(int d = 0; d < P_dim; d++) for(int dd = 0; dd < P_dim; dd++) P_season_terminal(t,s,a,d,dd) = P_t(d,dd);
      }
      P1 = P1 * get_P(n_regions, n_fleets, a, y, s, t, fleet_regions, cum_n_mu, n_mu, mu_row, mu_col, mu_pointer, mig_type(s), fracyr_seasons(t), FAA, MAA, mu, L);
      cum_n_mu += n_mu(s,y,t,a);
    }
    for(int i = 0; i < P_dim; i++) for(int j = 0; j < P_dim; j++) all_P(s,y,a,i,j) = P1(i,j);
  }
  //std::cout << "here 5"  << "\n";

  //get SSB from NAA_SSB, waa and mature
  SSB.setZero();
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++)
  {
    for(int a = 0; a < n_ages; a++) SSB(y,s) += NAA_SSB(s,y,a) * waa(waa_pointer_ssb(s)-1,y,a) * mature(s,y,a);
  }
  
  //get predicted numbers at age after year 1.
  pred_NAA.setZero();
  for(int s = 0; s < n_stocks; s++) 
  {
    for(int y = 1; y < n_years; y++)
    {
      for(int r = 0; r < n_regions; r++) 
      {
        if(r == spawn_regions(s)-1) //recruiting only in the stock region
        {
          if(recruit_model == 1) pred_NAA(s,y,0,r) = NAA(s,y-1,0,r); //random walkNAA(y,1)
          else
          {
            if(recruit_model == 2) pred_NAA(s,y,0,r) = exp(mean_rec_pars(s,0)); //random about mean
            else //BH stock recruit
            {
              if(recruit_model == 3) //BH stock recruit
              {
                pred_NAA(s,y,0,r) = exp(mean_rec_pars(s,0) + log(SSB(y-1,s)) - log(one + exp(mean_rec_pars(s,1)) * SSB(y-1,s)));
              }
              else //Ricker stock recruit
              {
                pred_NAA(s,y,0,r) = exp(mean_rec_pars(s,0) + log(SSB(y-1,s)) - exp(mean_rec_pars(s,1)) * SSB(y-1,s)); 
              }
            }
          }
        }
        else pred_NAA(s,y,0,r) = zero; // no recruitment for stock s in other regions.
      }
      //probs of  movement and survival and catch
      vector<Type> t_N(n_regions), t_pred_N(n_regions);
      for(int a = 1; a < n_ages; a++) 
      {
        for(int i = 0; i < n_regions; i++)
        {
          for(int j = 0; j < n_regions; j++)
          {
           pred_NAA(s,y,a,i) += all_P(s,y-1,a-1,j,i) * NAA(s,y-1,a-1,j);
          }
        }
        /*if(a == 1 && y == 4)
        {
          std::cout << "s: " << s << ", y: " << y << ", a: " << a << "\n";
          std::cout << NAA(s,y-1,a-1,0) << " " << NAA(s,y-1,a-1,1)  << "\n";
          std::cout << pred_NAA(s,y,a,0) << " " << pred_NAA(s,y,a,1)  << "\n";
          std::cout << all_P(s,y-1,a-1,0,0) << " " << all_P(s,y-1,a-1,0,1)  << "\n";
          std::cout << all_P(s,y-1,a-1,1,0) << " " << all_P(s,y-1,a-1,1,1)  << "\n";
        }*/
      }
      //add in extra for plus group
      for(int i = 0; i < n_regions; i++) for(int j = 0; j < n_regions; j++) 
      {
        pred_NAA(s,y,n_ages-1,i) += all_P(s,y-1,n_ages-1,j,i) * NAA(s,y-1,n_ages-1,j);
      }
    }
  }

  //need to be careful here about log_NAA random effects in regions where there will be zero predicted due to migration parameterization
  Type nll_NAA = zero;
  for(int s = 0; s < n_stocks; s++) for(int y = 1; y < n_years; y++) for(int a = 0; a < n_ages; a++) for(int r = 0; r < n_regions; r++) 
  {
    if(NAA_re_indicator(s,y-1,a,r) == 1)
    { 
     /*if(a == 1)
      {
        std::cout << "s: " << s << ", y: " << y << ", a: " << a << ", r: " << r << "\n";
        std::cout << NAA(s,y,a,r) << " " << pred_NAA(s,y,a,r) << " " << dnorm(log(NAA(s,y,a,r)),log(pred_NAA(s,y,a,r)),exp(log_NAA_sigma(NAA_sigma_pointers(s,a,r)-1)),1) << "\n";
      }*/
      nll_NAA -= dnorm(log(NAA(s,y,a,r)),log(pred_NAA(s,y,a,r)),exp(log_NAA_sigma(NAA_sigma_pointers(s,a,r)-1)),1);
    }
  }
  std::cout << "nll_NAA: " << nll_NAA << "\n";
  nll += nll_NAA;

  vector<Type> nll_agg_catch(n_fleets), nll_catch_acomp(n_fleets);
  nll_agg_catch.setZero();
  nll_catch_acomp.setZero();
  pred_CAA.setZero();
  pred_stock_CAA.setZero(); //(n_years,n_stocks,n_ages,n_regions);
  pred_stock_catch.setZero(); //(n_years,n_stocks,n_regions);
  pred_stock_prop_catch.setZero(); //(n_years,n_stocks,n_regions);
  pred_catch_region.setZero(); //(n_years,n_regions);
  pred_catch.setZero();

  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      //pred_catch(y,f) = zero;
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++) 
        {
          for(int rr = 0; rr < n_regions; rr++) if(fleet_regions(f)-1 == rr)
          {
            pred_stock_CAA(y,s,a,rr) += NAA(s,y,a,r) * all_P(s,y,a,r,n_regions + f);
            pred_stock_catch(y,s,rr) += waa(waa_pointer_fleets(f)-1,y,a) * NAA(s,y,a,r) * all_P(s,y,a,r,n_regions + f);
            pred_catch_region(y,rr) += NAA(s,y,a,r) * all_P(s,y,a,r,n_regions + f) * waa(waa_pointer_fleets(f)-1,y,a);
          }
          pred_CAA(y,f,a) +=  NAA(s,y,a,r) * all_P(s,y,a,r,n_regions + f);
        }
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      nll_agg_catch(f) -= dnorm(log(agg_catch(y,f)),log(pred_catch(y,f)),agg_catch_sigma(y,f),1);
      if(any_fleet_age_comp(f) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_fleets(f));
        for(int j = 0; j < n_age_comp_pars_fleets(f); j++) 
        {
          acomp_pars(j) = catch_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_catch_paa(y,f) == 1) 
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
            t_pred_paa(a) = pred_catch_paa(y,f,a);
            t_paa(a) = catch_paa(f * n_years + y,a);
          }
          //std::cout << "for y = " << y << ", get_acomp_ll: " << get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f)) << "\n";
          nll_catch_acomp(f) -= get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
        }
      }
    }
  }
  for(int y = 0; y < n_years; y++) for(int s = 0; s < n_stocks; s++) for(int r = 0; r < n_regions; r++)
  {    
    pred_stock_prop_catch(y,s,r) = pred_stock_catch(y,s,r)/pred_catch_region(y,r);
  }
  
  std::cout << "nll_agg_catch: " << "\n" << nll_agg_catch << "\n";
  nll += nll_agg_catch.sum();
  std::cout << "nll_catch_acomp: " << "\n" << nll_catch_acomp << "\n";
  nll += nll_catch_acomp.sum();

  vector<Type> nll_agg_indices(n_indices), nll_index_acomp(n_indices);
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  pred_indices.setZero();
  pred_IAA.setZero();
  for(int y = 0; y < n_years; y++) 
  {
    int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++) 
    {
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        //get numbers at age a for stock s in each region at time of spawning
        for(int s = 0; s < n_stocks; s++) pred_IAA(y,i,a) += QAA(y,i,a) * NAA_index(s,i,y,a);
        if(units_indices(i) == 1) pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);    
      }
      for(int a = 0; a < n_ages; a++) 
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }
      if(use_indices(y,i) == 1)
      {
        nll_agg_indices(i) -= dnorm(log(agg_indices(y,i)),log(pred_indices(y,i)),agg_index_sigma(y,i),1);
      }
      if(any_index_age_comp(i) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_indices(i));
        for(int j = 0; j < n_age_comp_pars_indices(i); j++) 
        {
          acomp_pars(j) = index_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_index_paa(y,i) > 0)
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
            t_pred_paa(a) = pred_index_paa(y,i,a);
            t_paa(a) = index_paa(i * n_years + y,a);
          }
          /*if(i == 1) 
          {
            Type acomp_ll = get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
            std::cout << "y: " << y << "\n";
            std::cout << index_Neff(y,i) << "\n";
            std::cout << age_comp_model_indices(i) << "\n";
            std::cout << "acomp_pars: " << "\n" << acomp_pars << "\n";
            std::cout << "index_aref(y,i): " << index_aref(y,i) << "\n";
            std::cout << "t_paa: " << "\n" << t_paa << "\n";
            std::cout << "t_pred_paa: " << "\n" << t_pred_paa << "\n";
            std::cout << "acomp_ll: " << acomp_ll << "\n";
          }*/
          nll_index_acomp(i) -= get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
        }
      }
    }
  }
  /*
  std::cout << "QAA: " << "\n";
  for(int i = 0; i < n_indices; i++)
  {
  std::cout << "index: " << i << "\n";
    for(int y = 0; y < n_years; y++)
    {
      for(int a = 0; a < n_ages; a++)
      {
        std::cout << QAA(y,i,a) << " ";
      }
      std::cout << "\n";
    }
  }
  std::cout << "pred_IAA: " << "\n";
  for(int i = 0; i < n_indices; i++)
  {
  std::cout << "index: " << i << "\n";
    for(int y = 0; y < n_years; y++)
    {
      for(int a = 0; a < n_ages; a++)
      {
        std::cout << pred_IAA(y,i,a) << " ";
      }
      std::cout << "\n";
    }
  }
  //std::cout << "NAA_index: " << "\n";
  //std::cout << NAA_index << "\n";
  */
  /*std::cout << "pred_indices: " << "\n";
  std::cout << pred_indices << "\n";
  */
  std::cout << "nll_agg_indices: " << "\n" << nll_agg_indices << "\n";
  nll += nll_agg_indices.sum();
  std::cout << "nll_index_acomp: " << "\n" << nll_index_acomp << "\n";
  nll += nll_index_acomp.sum();
  
  for(int s = 0; s < n_stocks; s++) for(int y = 0; y < n_years; y++) log_SSB(y,s) = log(SSB(y,s));
  
  //if(reportMode==0){
    REPORT(NAA);
    REPORT(pred_NAA);
    REPORT(SSB);
    REPORT(selblocks);
    REPORT(q);
    REPORT(F);
    REPORT(pred_catch);
    REPORT(pred_indices);
    REPORT(NAA_index);
    REPORT(QAA);
    ADREPORT(log_SSB);
    REPORT(pred_CAA);
    REPORT(pred_stock_CAA);
    REPORT(pred_stock_catch);
    REPORT(pred_stock_prop_catch);
    REPORT(pred_catch_region);
    REPORT(all_P);
    REPORT(P_season_terminal);
  //}
  
  return nll;
}

