template<class Type>
array<Type> get_FAA(matrix<Type>F, matrix<int> fleet_seasons, vector<matrix<Type>> selAA, matrix<int> selblock_pointer, int n_ages, int n_seasons){
  int n_fleets = F.dim(1);
  int n_y = F.dim(0);
  array<Type> FAA(n_fleets,n_y,n_seasons,n_ages);
  FAA.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) 
  {
    for(int t = 0; t < n_seasons; t++) if(fleet_seasons(f,t) == 1) 
    {
      FAA(f,y,t,a) = F(y,f) * selAA(selblock_pointer(y,f)-1)(y,a);
    }
  }
  return(FAA);
}
//done

template<class Type>
array<Type> get_FAA_tot(array<Type> FAA, vector<int> fleet_regions, int n_regions){
  int n_fleets = FAA.dim(0);
  int n_seasons = FAA.dim(2);
  int n_ages = FAA.dim(3);
  int n_y = FAA.dim(1);
  array<Type> FAA_tot(n_regions, n_y, n_seasons, n_ages);
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++) 
  {
    for(int t = 0; t < n_seasons; t++) if(fleet_seasons(f,t) == 1) 
    {
      FAA_tot(fleet_regions(f)-1,y,t,a) += FAA(f,y,t,a);
    }
  }
  return(FAA_tot);
}
//done

template <class Type>
array<Type> get_sel_proj(int y, array<Type> FAA, vector<int> avg_years_ind,
  vector<int> which_F_fleet, vector<int> which_F_season, vector<int> which_F_age){
    /* 
     get selectivity to project for next time step
                   y:  year of projection (>n_years_model)
                 FAA:  FAA array from main code.
       avg_years_ind:  integer vector of years to average F for projection
       which_F_fleet:  define which fleet has max F
      which_F_season:  define which season has max F
         which_F_age:  define which age has max F
    */
  //average F by fleet, season, and age is used to find selectivity (fleet,season,age) to project 
  //full F is the FAA for fleet, season and age defined by which_F_fleet,which_F_season, which_F_age
  int n_toavg = avg_years_ind.size();
  int n_fleets = FAA.dim(0);
  int n_seasons = FAA.dim(2);
  int n_ages = FAA.dim(3);

  array<Type> FAA_avg(n_fleets, n_seasons, n_ages);
  FAA_avg.setZero();
  for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++)
  {
    for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
      FAA_avg(f,t,a) += FAA(f,avg_years_ind(i),t,a)/Type(n_toavg);
    }
  }
  //get selectivity using average over avg.yrs
  array<Type> sel_proj(n_fleets,n_seasons,n_ages);
  //fully selected F across regions, seasons, and ages
  Type F_full = FAA_avg(which_F_fleet(y)-1,which_F_season(y)-1,which_F_age(y)-1);
  for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++){
    for(int a = 0; a < n_ages; a++) {
      sel_proj(f,t,a) = FAA_avg(f,t,a)/F_full;
    }
  }
  return(sel_proj);
}
//done

template <class Type>
array<Type> get_FAA_proj(int y, vector<int> proj_F_opt, array<Type> sel_proj, 
  array<Type> FAA, array<Type> NAA, array<Type> MAA, 
  array<Type> mature, array<Type> waa, vector<int> waa_pointer_totcatch, vector<int> waa_pointer_ssb, 
  matrix<Type> fracyr_SSB, matrix<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, vector<int> which_F_fleet, vector<int> which_F_season, vector<int> which_F_age, 
  Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, Type F_init, 
  matrix<Type> log_a, matrix<Type> log_b, vector<int> recruit_model, Type percentFMSY){
    /* 
     get FAA to project for next time step
                   y:  year of projection (>n_years_model)
          proj_F_opt:  for each year, how to specify F for projection (1 to 6)
                 FAA:  FAA array from main code.
                 NAA:  NAA array from main code
                 MAA:  MAA array from main code.
              mature:  maturity array from main code.
                 waa:  weight at age array 
waa_pointer_totcatch:  (n_regions) pointer for waa to use for tot catch (use function get_waa_y defined in helper.hpp)
     waa_pointer_ssb:  (n_stocks) pointer for waa to use for SSB (use function get_waa_y defined in helper.hpp)
          fracyr_SSB:  (n_stocks x n_years) vector of yearly fractions of the year when spawning occurs
            log_SPR0:  matrix (n_stocks x n_years) of yearly log(unfished SSB/R) 
       avg_years_ind:  integer vector of years to average F for projection
       n_years_model:  number of years before projection begins
       which_F_fleet:  define which fleet has max F
      which_F_season:  define which season has max F
         which_F_age:  define which age has max F
          percentSPR:  percentage (0-100) of unfished spawning potential to determine F_percentSPR
         proj_Fcatch:  vector (n_years_proj) of user specified Fishing mortality rates to project
        percentFXSPR:  percentage (0-100) of F_percentSPR to use in catch, e.g. GOM cod uses F = 75% F_40%SPR
              F_init:  initial value to use for FXSPR or FMSY newton method
               log_a:  (n_stocks x n_years) annual log(a) for stock-recruit relationship
               log_b:  (n_stocks x n_years) annual log(b) for stock-recruit relationship
       recruit_model:  (n_stocks) integer for which type of recruit model is assumed (= 3 or 4 for using Fmsy)
         percentFMSY:  percentage (0-100) of FMSY to use in catch.
    */
  int n_fleets = FAA.dim(0);
  int n_seasons = FAA.dim(2);
  int n_ages = FAA.dim(3);
  matrix<Type> waacatch = get_waa_y(waa, y, n_ages, waa_pointer_totcatch);
  matrix<Type> waassb = get_waa_y(waa, y, n_ages, waa_pointer_ssb);

  int n_toavg = avg_years_ind.size();
  int n_ages = waacatch.size();
  int proj_F_opt_y = proj_F_opt(y-n_years_model);

  //proj_F_opt == 1, last year F (default)
  array<Type> FAA_proj(n_fleets, n_seasons, n_ages);
  array<Type> FAA_tot_proj(n_regions,n_seasons, n_ages);
  FAA_tot_proj.setZero();
  FAA_proj.setZero();
  if(proj_F_opt_y == 1){ // last year F (default)
    for(int f = 0; f < n_fleets; f++) for(int t = 0; t < n_seasons; t++) for(int a = 0; a < n_ages; a++) {
      FAA_proj(f,t,a) = FAA(f,n_years_model-1,t,a); 
    }
  }
  else { //proj_F_opt_y>1
    //option 2: average F is by fleet and Ftot is sum of fleet averages
    //when there is more than 1 fleet, the sum of predicted catch across fleets will not generally equal the total catch using FAA_tot and waa_totcatch.
    Type F_full_proj = 0.0;
    if(proj_F_opt_y == 2) { //F_full is the same as that used to generate selectivity to project
      F_full_proj = FAA(which_F_fleet(y)-1, which_F_season(y)-1,which_F_age(y)-1); 
    }

    if(proj_F_opt_y == 4){ // user-specified F
      /*if(proj_Fcatch(y-n_years_model) < 1e-10){ // if F = 0, sel_proj is NaN
        FAA_proj.setZero();
      } else { */
        F_full_proj = Type(proj_Fcatch(y-n_years_model));
      //}
    }
    
    array<Type> sel_proj = get_sel_proj(y, FAA, avg_years_ind, which_F_fleet, which_F_season, which_F_age);
    /*
    if(proj_F_opt_y == 3){ // F at X% SPR
      F_full_proj = get_FXSPR(MAA_y, sel_proj, waassb, mat_y, percentSPR, fracyr_SSB_y, log_SPR0_y, F_init) * 0.01* percentFXSPR;
    }
    if(proj_F_opt_y == 5){ // calculate F from user-specified catch
      Type thecatch = proj_Fcatch(y-n_years_model);
      if(thecatch < 1e-10){ // if catch = 0, F = 0 and sel_proj is NaN
        //F_full_proj = 0.0; //already done
        //FAA_proj.setZero();
      } else {
        F_full_proj = get_F_from_log_Catch(thecatch, NAA_y, MAA_y, sel_proj, waacatch, F_init);
      }
    }
    /*if(proj_F_opt_y == 6){ //Fmsy
      vector<Type> log_a_y = log_a.row(y);
      vector<Type> log_b_y = log_b.row(y);
      F_full_proj = get_FMSY(log_a_y, log_b_y, MAA_y, sel_proj, waacatch, waassb, mat_y, fracyr_SSB_y, log_SPR0_y, recruit_model, F_init) * 0.01* percentFMSY;
    } */
    FAA_proj = F_full_proj * sel_proj;
  }
  return(FAA_proj);
}
