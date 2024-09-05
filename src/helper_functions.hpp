#define see(object) std::cout << #object ":\n" << object << "\n";

template <class Type>
Type square(Type x){return x*x;}

// general inverse logit fn
// template <class Type>
// Type geninvlogit(Type x, Type lo = 0.0, Type hi = 1.0, Type scale = 1.0){return lo + (hi-lo) * invlogit( scale * x);}

// general inverse logit fn
template <class T>
T geninvlogit(T x, T lo = 0.0, T hi = 1.0, T scale = 1.0){return lo + (hi-lo) * invlogit( scale * x);}

//VECTORIZE4_tttt(invlogit)

template <class Type>
Type logit(Type y, Type lo = 0.0, Type hi = 1.0, Type scale = 1.0){return (log(y - lo) - log(hi - y))/scale;}

//VECTORIZE4_tttt(logit);

// transformation to ensure correlation parameters are between -1 and 1
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}
VECTORIZE1_t(rho_trans)

template<class Type>
matrix<Type> extract_matrix_array3(array<Type> a, int index){ //matrix has to be the last two dims of the 3d array
  int nr = a.dim(1);
  int nc = a.dim(2);
  matrix<Type> mat(nr, nc);
  for(int i = 0; i < nr; i++) for(int j = 0; j < nc; j++) mat(i,j) = a(index,i,j);
  return mat;
}

template <class Type>
matrix<Type> get_avg_waa(array<Type> waa, vector<int> years, vector<int> pointers){
  
  int n_ages = waa.dim(2);
  int n_waa = pointers.size();
  int n_y = years.size();
  matrix<Type> avg_waa(n_waa,n_ages);
  avg_waa.setZero();
  for(int i = 0; i < n_waa; i++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++){
    avg_waa(i,a) += waa(pointers(i)-1, years(y), a)/Type(n_y);
  } 
  return avg_waa;
}

template <class Type>
array<Type> get_avg_waa_as_array(array<Type> waa, vector<int> years, vector<int> pointers){
  
  int n_ages = waa.dim(2);
  int n_waa = pointers.size();
  int n_y = years.size();
  array<Type> avg_waa(n_waa,n_ages);
  avg_waa.setZero();
  for(int i = 0; i < n_waa; i++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++){
    avg_waa(i,a) += waa(pointers(i)-1, years(y), a)/Type(n_y);
  } 
  return avg_waa;
}

template <class Type>
matrix<Type> get_avg_mat(array<Type> mat, vector<int> years){
  
  int n_ages = mat.dim(2);
  int n_stocks = mat.dim(0);
  int n_y = years.size();
  matrix<Type> avg_mat(n_stocks,n_ages);
  avg_mat.setZero();
  for(int i = 0; i < n_stocks; i++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++){
    avg_mat(i,a) += mat(i, years(y), a)/Type(n_y);
  } 
  return avg_mat;
}

template <class Type>
array<Type> get_avg_mat_as_array(array<Type> mat, vector<int> years){
  
  int n_ages = mat.dim(2);
  int n_stocks = mat.dim(0);
  int n_y = years.size();
  array<Type> avg_mat(n_stocks,n_ages);
  avg_mat.setZero();
  for(int i = 0; i < n_stocks; i++) for(int y = 0; y < n_y; y++) for(int a = 0; a < n_ages; a++){
    avg_mat(i,a) += mat(i, years(y), a)/Type(n_y);
  } 
  return avg_mat;
}

//extract a weight at age or maturity matrix for year y
template <class Type>
matrix<Type> get_matrix_y(array<Type> mat, int year){
  
  int n_ages = mat.dim(2);
  int n_mat = mat.dim(0); //n_stocks, n_fleets,...
  matrix<Type> avg_mat(n_mat,n_ages);
  avg_mat.setZero();
  for(int i = 0; i < n_mat; i++) for(int a = 0; a < n_ages; a++){
    avg_mat(i,a) += mat(i, year, a);
  } 
  return avg_mat;
}

//extract a weight at age or maturity matrix(array) for year y
template <class Type>
array<Type> get_matrix_y(array<Type> mat, int year){
  
  int n_ages = mat.dim(2);
  int n_mat = mat.dim(0); //n_stocks, n_fleets,...
  array<Type> avg_mat(n_mat,n_ages);
  avg_mat.setZero();
  for(int i = 0; i < n_mat; i++) for(int a = 0; a < n_ages; a++){
    avg_mat(i,a) += mat(i, year, a);
  } 
  return avg_mat;
}

template<class Type>
vector<Type> get_avg_ssbfrac(matrix<Type> ssbfrac, vector<int> years) {
  vector<Type> avg_ssbfrac(ssbfrac.cols());
  avg_ssbfrac.setZero();
  for(int s = 0; s < ssbfrac.cols(); s++) for(int y = 0; y < years.size(); y++) avg_ssbfrac(s) += ssbfrac(years(y),s)/Type(years.size());
  return avg_ssbfrac;
}

template<class Type>
matrix<Type> get_avg_SR_ab(matrix<Type> log_a, matrix<Type> log_b, vector<int> years_SR_ab, int do_log) {
  matrix<Type> SR_ab_avg(log_a.cols(), 2); // n_stocks x 2
  SR_ab_avg.setZero();
  //get average inputs over specified years
  for(int s = 0; s < log_a.cols(); s++) {
    for(int y = 0; y < years_SR_ab.size(); y++) {
      SR_ab_avg(s,0) += exp(log_a(y,s))/Type(years_SR_ab.size());
      SR_ab_avg(s,1) += exp(log_b(y,s))/Type(years_SR_ab.size());
    }
    if(do_log) {
      SR_ab_avg(s,0) = log(SR_ab_avg(s,0));
      SR_ab_avg(s,1) = log(SR_ab_avg(s,1));
    }
  }
  
  return(SR_ab_avg); //n_stocks x 2)
}

template <class Type>
matrix<Type> get_avg_FAA(array<Type> FAA, vector<int> years, int do_log){
  
  int n_fleets = FAA.dim(0);
  int n_ages = FAA.dim(2);
  matrix<Type> FAA_avg(n_fleets, n_ages);

  FAA_avg.setZero();
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){
    for(int y = 0; y < years.size(); y++) FAA_avg(f,a) += FAA(f,years(y),a)/Type(years.size()); //average F at fleet,season,age over years
    if(do_log) FAA_avg(f,a) = log(FAA_avg(f,a));
  }
  return FAA_avg;
}

template <class Type>
array<Type> get_avg_FAA_as_array(array<Type> FAA, vector<int> years, int do_log){
  
  int n_fleets = FAA.dim(0);
  int n_ages = FAA.dim(2);
  array<Type> FAA_avg(n_fleets, n_ages);

  FAA_avg.setZero();
  for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++){
    for(int y = 0; y < years.size(); y++) FAA_avg(f,a) += FAA(f,years(y),a)/Type(years.size()); //average F at fleet,season,age over years
    if(do_log) FAA_avg(f,a) = log(FAA_avg(f,a));
  }
  return FAA_avg;
}

template <class Type>
matrix<Type> get_avg_fleet_sel(array<Type> FAA, vector<int> avg_years_ind,
  int which_F_age){
    /* 
     get average selectivity. Typically to define referene points or for projections
                 FAA:  FAA (n_fleets x n_years x n_ages) array from main code.
       avg_years_ind:  integer vector of years to average FAA
         which_F_age:  define which age has max F
    */
  //average F by fleet, and age is used to find selectivity (fleet,season,age) to project 
  //full F is the FAA for fleet, season and age defined by which_F_age
  int n_toavg = avg_years_ind.size();
  int n_fleets = FAA.dim(0);
  int n_ages = FAA.dim(2);
  matrix<Type> FAA_avg(n_fleets, n_ages);
  FAA_avg.setZero();
  for(int f = 0; f < n_fleets; f++) {
    for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
      FAA_avg(f,a) += FAA(f,avg_years_ind(i),a)/Type(n_toavg);
    }
  }
  vector<Type> FAA_avg_tot = FAA_avg.colwise().sum();
  Type F_full = FAA_avg_tot(which_F_age-1);

  //get selectivity using average over avg.yrs
  matrix<Type> sel(n_fleets,n_ages);
  //fully selected F across regions, seasons, and ages
  for(int f = 0; f < n_fleets; f++){
    for(int a = 0; a < n_ages; a++) {
      sel(f,a) = FAA_avg(f,a)/F_full;
    }
  }
  return sel;
}

template <class Type>
array<Type> get_avg_fleet_sel_as_array(array<Type> FAA, vector<int> avg_years_ind,
  int which_F_age, int by_fleet = 0){
    /* 
     get average selectivity. Typically to define referene points or for projections
                 FAA:  FAA (n_fleets x n_years x n_ages) array from main code.
       avg_years_ind:  integer vector of years to average FAA
         which_F_age:  define which age has max F
    */
  //average F by fleet, and age is used to find selectivity (fleet,season,age) to project 
  //full F is the FAA for fleet, season and age defined by which_F_age
  int n_toavg = avg_years_ind.size();
  int n_fleets = FAA.dim(0);
  int n_ages = FAA.dim(2);
  matrix<Type> FAA_avg(n_fleets, n_ages);
  FAA_avg.setZero();
  for(int f = 0; f < n_fleets; f++) {
    for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
      FAA_avg(f,a) += FAA(f,avg_years_ind(i),a)/Type(n_toavg);
    }
  }
  vector<Type> FAA_avg_tot = FAA_avg.colwise().sum();
  Type F_full = FAA_avg_tot(which_F_age-1);

  //get selectivity using average over avg.yrs
  array<Type> sel(n_fleets,n_ages);
  //fully selected F across regions, seasons, and ages
  for(int f = 0; f < n_fleets; f++){
    for(int a = 0; a < n_ages; a++) {
      if(by_fleet == 0) sel(f,a) = FAA_avg(f,a)/F_full;
      else sel(f,a) = FAA_avg(f,a)/ FAA_avg(f,which_F_age-1);
    }
  }
  return sel;
}

template <class Type>
array<Type> get_avg_M(array<Type> log_M, vector<int> years, int do_log){
  
  int n_ages = log_M.dim(3);
  int n_regions = log_M.dim(1);
  int n_stocks = log_M.dim(0);
  int n_y = years.size();
  array<Type> avg_M(n_stocks,n_regions, n_ages);
  avg_M.setZero();
  for(int i = 0; i < n_stocks; i++) for(int r = 0; r < n_regions; r++) for(int a = 0; a < n_ages; a++){
    for(int y = 0; y < n_y; y++) avg_M(i, r, a) += exp(log_M(i, r, years(y), a))/Type(n_y);
    if(do_log) avg_M(i, r, a) = log(avg_M(i, r, a));
  }
  return avg_M;
}


template <class Type>
vector<Type> get_avg_L(matrix<Type> L, vector<int> years, int do_log){
  int n_y = years.size();
  int n_regions = L.cols();
  vector<Type> avg_L(n_regions);
  avg_L.setZero();
  for(int i = 0; i < n_regions; i++) {
    for(int y = 0; y < n_y; y++) avg_L(i) += L(years(y),i)/Type(n_y);
    if(do_log) avg_L(i) = log(avg_L(i));
  }
  return avg_L;
}


template <class Type>
vector<Type> get_waa_y(array<Type> waa, int y, int na, int pointer){
  vector<Type> waay(na);
  for(int a = 0; a < na; a++) waay(a) = waa(pointer-1, y, a);
  return(waay);
}

template <class Type>
matrix<Type> get_waa_y(array<Type> waa, int y, int na, vector<int> pointer){
  matrix<Type> waay(pointer.size(),na);
  for(int i = 0; i < pointer.size(); i++) for(int a = 0; a < na; a++) waay(i,a) = waa(pointer(i)-1, y, a);
  return(waay);
}


template <class Type>
matrix<Type> get_waacatch_y(array<Type> waa, int y, int na, vector<int> pointer){
  matrix<Type> waay(pointer.size(), na);
  for(int f = 0; f< pointer.size(); f++) for(int a = 0; a < na; a++) waay(f,a) = waa(pointer(f)-1, y, a);
  return(waay);
}

// transform vector 'x' into matrix of orthogonal polynomials, with degree/ncols = 'degree'
// note that the # datapoints, length(x), must be greater than 'degree' - this is checked on the R side before calling TMB
// degree assumed to be > 1 (could call this for linear Ecov models but slows down code, so don't)
// algorithm from https://stackoverflow.com/questions/39031172/how-poly-generates-orthogonal-polynomials-how-to-understand-the-coefs-ret/39051154
template <class Type>
matrix<Type> poly_trans(vector<Type> x, int degree, int n_years_model, int n_years_proj)
{
  vector<Type> x_model = x.head(n_years_model);
  Type x_mean = x_model.mean();
  vector<Type> x_centered = x_model - x_mean;
  
  vector<Type> beta(degree);
  vector<Type> alpha(degree);
  vector<Type> norm2(degree);
  beta.setZero();
  alpha.setZero();
  norm2.setZero();
  matrix<Type> X(n_years_model, degree);
  X = x_centered.replicate(1,degree);
  
  Type new_norm = (x_centered * x_centered).sum();
  norm2(0) = new_norm;
  alpha(0) = (x_centered * x_centered * x_centered).sum() / new_norm;
  beta(0) = new_norm / n_years_model;
  
  // degree 2 (assume degree > 1)
  Type old_norm = new_norm;
  vector<Type> Xi(n_years_model);
  Xi = (x_centered - alpha(0)).array() * X.col(0).array() - beta(0);
  X.col(1) = Xi;
  vector<Type> tmp2 = Xi * Xi;
  new_norm = tmp2.sum();
  norm2(1) = new_norm;
  alpha(1) = (tmp2 * x_centered).sum() / new_norm;
  beta(1) = new_norm / old_norm;
  old_norm = new_norm;
  
  // for degrees > 2
  if(degree > 2){
    for(int i=3; i<degree+1; i++){
      Xi =  (x_centered - alpha(i-2)).array() * X.col(i-2).array() - beta(i-2)*X.col(i-3).array();
      X.col(i-1) = Xi;
      new_norm = (Xi * Xi).sum();
      norm2(i-1) = new_norm;
      alpha(i-1) = (Xi * Xi * x_centered).sum() / new_norm;
      beta(i-1) = new_norm / old_norm;
      old_norm = new_norm;      
    }
  }

  // scale
  vector<Type> scale = sqrt(norm2);
  matrix<Type> finalX(n_years_model + n_years_proj, degree);
  for(int j=0; j<degree; j++){
    for(int i=0; i<n_years_model; i++){
      finalX(i,j) = X(i,j) / scale(j);
    }
  }

  // do proj years if necessary
  if(n_years_proj > 0){
    vector<Type> x_proj = x.tail(n_years_proj);
    vector<Type> x_centered_proj = x_proj - x_mean;
    matrix<Type> X_proj(n_years_proj, degree);
    X_proj = x_centered_proj.replicate(1,degree);

    // degree 2
    vector<Type> Xi_proj(n_years_proj);
    Xi_proj = (x_centered_proj - alpha(0)).array() * X_proj.col(0).array() - beta(0);
    X_proj.col(1) = Xi;

    // degree > 2
    if(degree > 2){
      for(int i=3; i<degree+1; i++){
        Xi_proj =  (x_centered_proj - alpha(i-2)).array() * X_proj.col(i-2).array() - beta(i-2)*X_proj.col(i-3).array();
        X_proj.col(i-1) = Xi;      
      }
    }

    for(int j=0; j<degree; j++){
      for(int i=0; i<n_years_proj; i++){
        finalX(n_years_model+i,j) = X_proj(i,j) / scale(j);
      }
    }
  }
  
  return finalX;
}

template <class Type>
Type d2dar1(array<Type> delta, Type tf_rho_r, Type tf_rho_c, vector<Type> log_sig_c, int use_dns=1){
  
  Type rho_r = geninvlogit(tf_rho_r, Type(-1), Type(1), Type(1)); //rho_year
  Type rho_c = geninvlogit(tf_rho_c,Type(-1), Type(1), Type(1)); //rho_age
  vector<Type> sig_c = exp(log_sig_c); //sd at age
  vector<Type> marg_sig = sig_c * pow((1 - pow(rho_r,2))*(1 - pow(rho_c,2)),-0.5); // marginal at age across years and ages
  vector<Type> sig_r1 = sig_c * pow(1-pow(rho_r,2),-0.5); //marginal sd at age across years
  Type sig_c1 = sig_c(0) * pow(1-pow(rho_c,2),-0.5); //marginal sd age 1 across ages given year
  Type res = 0;
  if(use_dns == 0){
    //first row
    res -= dnorm(delta(0,0), Type(0), marg_sig(0),1); //marginal across year and age
    for(int a = 1; a < delta.dim(1); a++){
      res -= dnorm(delta(0,a), rho_c * delta(0,a-1) * sig_c(a) / sig_c(a-1), sig_r1(a), 1); //marginal across years conditional on age a-1
    }
    //subsequent rows
    for(int y = 1; y < delta.dim(0); y++){
      res -= dnorm(delta(y,0), rho_r * delta(y-1,0), sig_c1, 1); //marginal at age 1 across ages and conditional on year y-1
    }
    for(int y = 1; y < delta.dim(0); y++) {
      for(int a = 1; a < delta.dim(1); a++) {
        res -= dnorm(delta(y,a), rho_r *delta(y-1,a) + rho_c * (delta(y,a-1) - rho_r * delta(y-1,a-1)) * sig_c(a) /sig_c(a-1), sig_c(a), 1); 
      }
    }
  } else{
    using namespace density;
    res = SEPARABLE(VECSCALE(AR1(rho_c), marg_sig), AR1(rho_r))(delta);
  }
  return res;
}

template <class Type>
array<Type> r2dar1(array<Type> delta, Type tf_rho_r, Type tf_rho_c, vector<Type> log_sig_c, int use_dns=1, int ystart=0, int bias_correct=0){
  
  Type rho_r = geninvlogit(tf_rho_r, Type(-1), Type(1), Type(1)); //rho_year (rows)
  Type rho_c = geninvlogit(tf_rho_c,Type(-1), Type(1), Type(1)); //rho_age (columns)
  vector<Type> sig_c = exp(log_sig_c); //conditional sd at age
  vector<Type> marg_sig = sig_c * pow((1 - pow(rho_r,2))*(1 - pow(rho_c,2)),-0.5); // marginal at age across years and ages
  vector<Type> sig_r1 = sig_c * pow(1-pow(rho_r,2),-0.5); //marginal sd at age across years
  Type sig_c1 = sig_c(0) * pow(1-pow(rho_c,2),-0.5); //marginal sd age 1 across ages given year
  // Type res = 0;
  array<Type> delta_out = delta;
  if(use_dns == 0){
    vector<Type> mu_c(marg_sig.size());
    mu_c.setZero();
    if(bias_correct) mu_c = - 0.5 * marg_sig * marg_sig;
    if(ystart == 0){ //first row, no conditioning on previous years
      delta_out(0,0) = rnorm(mu_c(0), marg_sig(0)); //marginal across year and age at age 1
      // if(bias_correct) delta_out(0,0) -= 0.5 * pow(marg_sig(0),2);
      for(int a = 1; a < delta.dim(1); a++){
        delta_out(0,a) = rnorm(mu_c(a) + (delta_out(0,a-1) - mu_c(a-1)) * rho_c * sig_c(a) / sig_c(a-1), sig_r1(a)); //marginal across years conditional on age a-1
        // if(bias_correct) delta_out(0,a) -= 0.5 * pow(marg_sig(a),2) * (1 - rho_c * sig_c(a-1)/sig_c(a));
      }
      for(int y = 1; y < delta.dim(0); y++){
        delta_out(y,0) = rnorm(rho_r * delta_out(y-1,0), sig_c1); //marginal at age 1 across ages and conditional on year y-1
        delta_out(y,0) += (1 - rho_r) * mu_c(0);
        // if(bias_correct) delta_out(y,0) -= 0.5 * pow(marg_sig(0),2) * (1-rho_r);
        for(int a = 1; a < delta.dim(1); a++) {
          delta_out(y,a) = rnorm(rho_r *delta_out(y-1,a) + rho_c * (delta_out(y,a-1) - rho_r * delta_out(y-1,a-1)) * sig_c(a) /sig_c(a-1), sig_c(a)); 
          delta_out(y,a) += (1 - rho_r) * (mu_c(a) - mu_c(a-1) * rho_c * sig_c(a)/sig_c(a-1));
          // if(bias_correct) delta_out(y,a) -= 0.5 * pow(marg_sig(a),2) * (1 - rho_r) * (1 - rho_c * sig_c(a-1)/sig_c(a));
        }
      }
    } else {
      //subsequent rows
      for(int y = ystart; y < delta.dim(0); y++){
        delta_out(y,0) = rnorm(rho_r * delta_out(y-1,0), sig_c1); //marginal at age 1 across ages and conditional on year y-1
        delta_out(y,0) += (1 - rho_r) * mu_c(0);
        // if(bias_correct) delta_out(y,0) -= 0.5 * pow(marg_sig(0),2) * (1-rho_r);
        for(int a = 1; a < delta.dim(1); a++) {
          delta_out(y,a) = rnorm(rho_r *delta_out(y-1,a) + rho_c * (delta_out(y,a-1) - rho_r * delta_out(y-1,a-1)) * sig_c(a) /sig_c(a-1), sig_c(a)); 
          delta_out(y,a) += (1 - rho_r) * (mu_c(a) - mu_c(a-1) * rho_c * sig_c(a)/sig_c(a-1));
          // if(bias_correct) delta_out(y,a) -= 0.5 * pow(marg_sig(a),2) * (1 - rho_r) * (1 - rho_c * sig_c(a-1)/sig_c(a));
        }
      }
    }
  } else{
    using namespace density;
    SEPARABLE(VECSCALE(AR1(rho_c), marg_sig),AR1(rho_r)).simulate(delta_out); // scaled here
    if(bias_correct) for(int y = 0; y < delta.dim(0); y++) for(int a = 0; a < delta.dim(1); a++) delta_out(y,a) -= 0.5 * pow(marg_sig(a),2);
  }
  return delta_out;
}

template <class Type>
Type dar1(vector<Type> delta, Type tf_rho, Type log_sig, int use_dns){
  
  Type rho = geninvlogit(tf_rho, Type(-1), Type(1), Type(1)); 
  Type sig = exp(log_sig); 
  Type marg_sig = sig * pow(1 - pow(rho,2),-0.5);
  Type res = 0;
  if(use_dns == 0){
    res -= dnorm(delta(0), Type(0), marg_sig,1); //marginal across year and age
    for(int y = 1; y < delta.size(); y++){
      res -= dnorm(delta(y), rho * delta(y-1), sig, 1); //marginal at age 1 and conditional on year y-1
    }
  } else{
    using namespace density;
    res = SCALE(AR1(rho), marg_sig)(delta);
  }
  return res;
}


template <class Type>
vector<Type> rar1(vector<Type> delta, Type tf_rho, Type log_sig, int use_dns, int ystart=0, int bias_correct=0){
  
  Type rho = geninvlogit(tf_rho, Type(-1), Type(1), Type(1)); 
  Type sig = exp(log_sig); 
  Type marg_sig = sig * pow(1 - pow(rho,2),-0.5);
  vector<Type> delta_out = delta;
  if(use_dns == 0){
    if(ystart == 0) {
      delta_out(0) = rnorm(Type(0), marg_sig); //marginal across year and age
      if(bias_correct) delta_out(0) -= 0.5*pow(marg_sig,2);
      for(int y = 1; y < delta.size(); y++){
        delta_out(y) = rnorm(rho * delta_out(y-1), sig); //marginal at age 1 and conditional on year y-1
        if(bias_correct) delta_out(y) -= 0.5*pow(marg_sig,2) * (1-rho);
      }
    } else {
      for(int y = ystart; y < delta.size(); y++) {
        delta_out(y) = rnorm(rho * delta_out(y-1), sig);
        if(bias_correct) delta_out(y) -= 0.5*pow(marg_sig,2) * (1-rho);
      }
    }
  } else {
    using namespace density;
    SCALE(AR1(rho), marg_sig).simulate(delta_out);
    if(bias_correct) delta_out -= 0.5*pow(marg_sig,2);
  }
  return delta_out;
}
