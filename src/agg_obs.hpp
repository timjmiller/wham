// Simulate from tweedie
// Adapted from tweedie::rtweedie function in R
template<class Type>
Type rTweedie( Type mu, Type phi, Type power){
 Type lambda = pow(mu, Type(2.0) - power) / (phi * (Type(2.0) - power));
 Type alpha = (Type(2.0) - power) / (Type(1.0) - power);
 Type gam = phi * (power - Type(1.0)) * pow(mu, power - Type(1.0));
 Type N = rpois(lambda);
 Type B = rgamma(N * alpha, gam);   /// Using Shape-Scale parameterization
 return B;
}

template<class Type>
struct gengamma_log_b_objective {
  /* Data and parameter objects for calculation: */
  Type k;
  Type CV;
  int trace=0;

  /* Constructor */
  gengamma_log_b_objective(
  Type k_,
  Type CV_,
  int trace_) :
    k(k_), CV(CV_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> log_b) { //find gengamma p (=beta in other papers) consistent with CV+, mu, lambda.
    T b = exp(log_b(0));
    vector<T> diff(1), diff2(1);
    diff(0) = (T(CV)*T(CV)+1.0) - exp(lgamma(T(k)) +lgamma(T(k)+2.0/b) - 2.0*lgamma(T(k) + 1.0/b));
    diff2(0) = diff(0) * diff(0);
    if(trace) {
      see(T(CV)*T(CV)+1.0);
      see(exp(lgamma(T(k)) +lgamma(T(k)+2.0/b) - 2.0*lgamma(T(k) + 1.0/b)));
      see(diff);
      see(b);
    }
    return diff2;
  }
};

template <class Type>
matrix<Type> get_gengamma_log_b(Type CV, Type k, int n = 30, int trace = 0){
  matrix<Type> log_b_iter(n,3);
  log_b_iter.setZero();
  gengamma_log_b_objective<Type> log_b_obj(k, CV, trace);
  vector<Type> log_b_i(1);
  log_b_i(0) = 0.0;
  for (int i=0; i<n-1; i++) {
    log_b_i(0) = log_b_iter(i,0); //current log_b
    log_b_iter(i,1) = log_b_obj(log_b_i)(0); //f(log_b)
    log_b_iter(i,2) = autodiff::jacobian(log_b_obj,log_b_i)(0,0); //f'(log_b)
    log_b_iter(i+1,0) = log_b_iter(i,0) - log_b_iter(i,1)/log_b_iter(i,2); //next log_b
  }
  return(log_b_iter);
}

template <class Type>
Type dgengamma(Type x, Type a, Type b, Type Q, int give_log=0, int prentice = 1){
  //a = scale
  //b = shape
  Type k = 1/(Q*Q); //Q = 1/sqrt(k), x goes to lognormal as Q goes to 0.
  Type logres = 0.0;
  if(!prentice) {
    logres = log(b) + (k*b-1) * log(x) - exp((log(x)-log(a))*b) - k*b*log(a) - lgamma(k); 
  } else {
    Type sigma = 1/(b*Q);
    Type mu = log(a) + log(k)/b;
    Type w = (log(x) - mu)/sigma;
    logres = -log(sigma) - log(x) + log(Q) *(1.0 - 2.0*k) + k * (Q*w - exp(Q*w)) - lgamma(k);
  }
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type pgengamma(Type x, Type a, Type b, Type Q, int give_log=0, int prentice = 1){
  //not concerned with Q<0
  Type k = 1/(Q*Q); //Q = 1/sqrt(k), x goes to lognormal as Q goes to 0.
  Type res = 0.0;
  if(!prentice) {
    Type y = pow(x/a,b);
    res = pgamma(y,k,Type(1.0));
  } else {
    Type sigma = 1/(b*Q);
    Type mu = log(a) + log(k)/b;
    Type w = (log(x) - mu)/sigma;
    Type expnu = exp(Q*w)*k; // = (x/a)^b
    res = pgamma(expnu,k,Type(1.0));
  }
  if(give_log) return log(res); else return res;
}

template<class Type>
Type rgengamma(Type a, Type b, Type Q){
  Type k = 1/(Q*Q); //Q = 1/sqrt(k), x goes to lognormal as Q goes to 0.
  Type w = log(rgamma(k,Type(1.0)));
  return a*exp(w/b);
}

template<class Type>
struct CPG_log_beta_objective {
  /* Data and parameter objects for calculation: */
  Type mu;
  Type lambda;
  Type cvplus;
  int trace;

  /* Constructor */
  CPG_log_beta_objective(
  Type mu_,
  Type lambda_,
  Type cvplus_,
  int trace_) :
    mu(mu_), lambda(lambda_), cvplus(cvplus_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> log_beta) { //find CPG beta consistent with CV+, mu, lambda.
    T beta = exp(log_beta(0));
    T alpha = T(mu) * beta /T(lambda);
    T p_0 = exp(-T(lambda));
    T power = (2.0 + alpha)/(1.0 + alpha);
    T phi = pow(mu,2.0-power) * (T(cvplus) * T(cvplus) + p_0);
    T beta_fn = pow(mu, 1.0 - power) / ((power - 1.0) * phi);
    vector<T> diff(1);
    diff(0) = log_beta(0) - log(beta_fn);
    if(trace) see(beta);    
    return diff * diff;
  }
};

template <class Type>
matrix<Type> get_CPG_log_beta(Type mu, Type lambda, Type cvplus, int n = 30, int trace = 0){
  matrix<Type> log_beta_iter(n,3);
  log_beta_iter.setZero();
  CPG_log_beta_objective<Type> beta_obj(mu, lambda, cvplus, trace);
  vector<Type> log_beta_i(1);
  for (int i=0; i<n-1; i++) {
    log_beta_i(0) = log_beta_iter(i,0);
    log_beta_iter(i,1) = beta_obj(log_beta_i)(0); //f(log_beta)
    log_beta_iter(i,2) = autodiff::jacobian(beta_obj,log_beta_i)(0,0); //f'(log_beta)
    ///Type jac_beta_obj = autodiff::jacobian(beta_obj,log_beta_i)(0,0);
    log_beta_iter(i+1,0) = log_beta_iter(i,0) - log_beta_iter(i,1)/log_beta_iter(i,2); 
    //log_beta_iter(i+1) = log_beta_iter(i) - beta_obj(log_beta_i)(0)/jac_beta_obj; 
  }
  return(log_beta_iter);
}

template <class Type>
matrix<Type> get_log_pred(matrix<Type> pred, matrix<Type> obs_pars, matrix<Type> obs_sd, vector<int> obs_model, int bias_correct_oe){
  int n_y = obs_sd.rows();
  int n_series = obs_sd.cols();
  matrix<Type> log_pred(n_y,n_series);
  log_pred.setZero();
  Type sd = 1.0;
  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++) {
    log_pred(y,i) = log(pred(y,i));
    if(bias_correct_oe){
      if((obs_model(i) == 1) | (obs_model(i) == 7)) {
        sd = obs_sd(y,i) * exp(obs_pars(i,0));
        log_pred(y,i) -= 0.5*exp(2*log(sd));
      }
      if((obs_model(i) == 2) | (obs_model(i) == 8)){
        sd = obs_sd(y,i) * (1.0 + exp(obs_pars(i,0)));
        log_pred(y,i) -= 0.5*exp(2*log(sd));
      }
    }
  }
  return log_pred;
}

template <class Type>
array<Type> get_annual_agg_obs_nll_pars(vector<int> agg_obs_model, matrix<Type> log_pred, matrix<Type> obs_sd, matrix<Type> effort, matrix<Type> agg_obs_pars, int n_iter = 30){
  int n_y = log_pred.rows();
  int n_series = log_pred.cols();
  array<Type> annual_nll_pars(n_series,n_y,4);
  annual_nll_pars.setZero();

  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++) {
    if((agg_obs_model(i) == 1) | (agg_obs_model(i) == 7)){ //normal for log-obs: up to 1 free parameters (scalar)  (delta: up to 3)
      annual_nll_pars(i,y,0) = log_pred(y,i); //mu for normal pdf of log-obs
      annual_nll_pars(i,y,1) = obs_sd(y,i)*exp(agg_obs_pars(i,0)); //sigma for normal pdf of log-obs
    }
    if((agg_obs_model(i) == 2) | (agg_obs_model(i) == 8)){ //normal for log-obs: with CV and SD(log(catch)) required to be greater than input  (delta: up to 3)
      annual_nll_pars(i,y,0) = log_pred(y,i); //mu for normal pdf of log-obs
      annual_nll_pars(i,y,1) = obs_sd(y,i) * (1.0 + exp(agg_obs_pars(i,0)));
    }
    if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 9)){ //gamma: up to 1 free parameters (scalar) (delta: up to 3)
      //gamma rate parameter is returned: beta = E(x)/(scalar * obs_sd)^2, (scalar = 1 means mu and CV determine distribution)
      // SD = E(x) / beta
      annual_nll_pars(i,y,1) = exp(log_pred(y,i) -2.0*(agg_obs_pars(i,0) + log(obs_sd(y,i)))); //gamma beta(rate) parameter
      annual_nll_pars(i,y,0) = annual_nll_pars(i,y,1) * exp(log_pred(y,i)); //gamma alpha (shape) parameter: beta * mu
    }
    if((agg_obs_model(i) == 4) | (agg_obs_model(i) == 10)){ //gamma:  with variance required to be > mu/beta  (delta: up to 3)
      //gamma rate parameter is returned: beta =  E(x)/((scalar+1) * obs_sd)^2,, CV and SD required to be greater than input
      // V = mu / beta
      annual_nll_pars(i,y,1) = exp(log_pred(y,i) -2.0*(log(exp(agg_obs_pars(i,0)) + 1.0) + log(obs_sd(y,i)))); //gamma beta(rate) parameter
      annual_nll_pars(i,y,0) = annual_nll_pars(i,y,1) * exp(log_pred(y,i)); //gamma alpha (shape) parameter: beta * mu
    }
    if((agg_obs_model(i) == 5) | (agg_obs_model(i) == 11)){ //generalized gamma: default 1 free parameters (delta: up to 3)
      //Type CV2 = exp(pow(obs_sd(y,i),2.0))-1.0; //CV^2
      //Type cv = pow(CV2,0.5);
      Type cv = obs_sd(y,i) * exp(-log_pred(y,i));
      Type k = exp(-2.0*agg_obs_pars(i,0)); // 1/(Q*Q), 
      matrix<Type> log_b_iter = get_gengamma_log_b(cv,k,n_iter); //iterations for finding log(b) of Generalized Gamma RVs.
      Type log_a = log_pred(y,i) +lgamma(k) - lgamma(k+exp(-log_b_iter(n_iter-1,0))); //theta = a in wikipedia
      annual_nll_pars(i,y,0) = exp(log_a); //scale par
      annual_nll_pars(i,y,1) = exp(log_b_iter(n_iter-1,0)); //shape par
      annual_nll_pars(i,y,2) = exp(agg_obs_pars(i,0));//Gamma par Q = 1/sqrt(k);
    }
    if(agg_obs_model(i) == 6){ //Tweedie (compound Poisson-Gamma): up to 2 free parameters (possibly some annual CV RE)
      Type CV2_pos = exp(pow(obs_sd(y,i),2.0))-1.0; //CV^2
      Type cv_pos = pow(CV2_pos,0.5);
      Type mu_pos = exp(log_pred(y,i));
      Type k1 = exp(agg_obs_pars(i,0)); //rate of positive obs with pred_obs (q*sum(sel*N)))
      Type k2 = exp(agg_obs_pars(i,1)); //rate of positive obs with measure of "effort")
      Type p_pos = (1 - exp(-k1*mu_pos)) * (1 - exp(-k2*effort(y,i))); //probability of positive obs: 1 - prob Poisson RV = 0.
      Type lambda = -log(1.0-p_pos); //mean of Poisson distributed number of positive gamma RVs.
      Type mu = mu_pos * p_pos; //mean of CPG
      matrix<Type> log_beta_iter = get_CPG_log_beta(mu, lambda, cv_pos); //iterations for finding log(beta) of Gamma RVs.
      Type alpha = mu * exp(log_beta_iter(log_beta_iter.rows()-1,0)) /lambda; //Gamma RV parameter
      Type power = (2.0 + alpha)/(1.0 + alpha); //Tweedie parameter: 1<power<2 for CPG
      Type phi = pow(mu,2.0-power) * (cv_pos*cv_pos + 1 - p_pos); //Tweedie parameter
      annual_nll_pars(i,y,0) = mu; //mean of Tweedie
      annual_nll_pars(i,y,1) = power;
      annual_nll_pars(i,y,2) = phi;
    }
    if(agg_obs_model(i) > 6){ //delta-models: up to 3 free parameter
      //7: delta-lognormal, 8: delta-lognormal(scalar>1), 9: delta-gamma, 10: delta-gamma(var>mu/beta), 11: delta-generalized gamma
      Type k1 = exp(agg_obs_pars(i,1)); //rate of positive obs with pred_obs (q*sum(sel*N)))
      Type k2 = exp(agg_obs_pars(i,2)); //rate of positive obs with measure of "effort")
      Type mu_pos = exp(log_pred(y,i));
      Type p_pos = (1 - exp(-k1*mu_pos)) * (1 - exp(-k2*effort(y,i))); //probability of positive obs: 1 - prob Poisson RV = 0.
      if(agg_obs_model(i) < 11) {
        annual_nll_pars(i,y,2) = p_pos; //probablity of positive obs
      } else{ //11: delta-generalized gamma
        annual_nll_pars(i,y,3) = p_pos; //probablity of positive obs
      }
    }
  }
  return annual_nll_pars;
}


template <class Type>
matrix<Type> get_nll_agg_obs(vector<int> agg_obs_model, matrix<int> obs_pos, array<Type> annual_nll_pars,
  vector<Type> obsvec, matrix<int> use_obs, matrix<int> keep_ind, data_indicator<vector<Type>, Type> keep){
  int n_series = annual_nll_pars.dim(0);
  int n_y = annual_nll_pars.dim(1);

  //for delta(hurdle) models and Tweedie, need to make sure call to TMB::oneStepPredict has discreteSupport = 0 and range = c(0,Inf)
  matrix<Type> nll_agg_obs(n_y,n_series);
  nll_agg_obs.setZero();

  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++) if(use_obs(y,i)){
    if((agg_obs_model(i) <3) | (agg_obs_model(i) == 7) | (agg_obs_model(i) == 8)) { //(delta) normal for log-obs, two types
      Type mu = annual_nll_pars(i,y,0);
      Type sig = annual_nll_pars(i,y,1);
      if(agg_obs_model(i) < 3){
        nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * dnorm(obsvec(keep_ind(y,i)), mu, sig,1);
        nll_agg_obs(y,i) -= keep.cdf_lower(keep_ind(y,i)) * log(squeeze(pnorm(obsvec(keep_ind(y,i)), mu, sig)));
        nll_agg_obs(y,i) -= keep.cdf_upper(keep_ind(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_ind(y,i)), mu, sig)));
      } else { //delta model
        Type p_pos = annual_nll_pars(i,y,2);
        if(obs_pos(y,i)) nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (log(p_pos) + dnorm(obsvec(keep_ind(y,i)), mu, sig,1));
        else nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * log(1-p_pos);
      }
    }
    if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 4) | (agg_obs_model(i) == 9) | (agg_obs_model(i) == 10))//(delta) gamma, two types
    if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 4)){ //not delta model
      Type alpha = annual_nll_pars(i,y,0);
      Type scale = 1/annual_nll_pars(i,y,1);
      if(agg_obs_model(i) < 5){
        nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * dgamma(obsvec(keep_ind(y,i)), alpha, scale,1);
        nll_agg_obs(y,i) -= keep.cdf_lower(keep_ind(y,i)) * log(squeeze(pgamma(obsvec(keep_ind(y,i)), alpha, scale)));
        nll_agg_obs(y,i) -= keep.cdf_upper(keep_ind(y,i)) * log(1.0 - squeeze(pgamma(obsvec(keep_ind(y,i)), alpha, scale)));
      } else { //delta model
        Type p_pos = annual_nll_pars(i,y,2);
        if(obs_pos(y,i)) nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (log(p_pos) + dgamma(obsvec(keep_ind(y,i)), alpha, scale,1));
        else nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * log(1-p_pos);
      }
    }
    if((agg_obs_model(i) == 5) | (agg_obs_model(i) == 11)){ //(delta) generalized gamma
      Type a = annual_nll_pars(i,y,0); //scale par
      Type b = annual_nll_pars(i,y,1); //shape par
      Type Q = annual_nll_pars(i,y,2); //Gamma par Q = 1/sqrt(k);
      if(agg_obs_model(i) == 5){
        nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * dgengamma(obsvec(keep_ind(y,i)), a, b, Q, 1, 1); //Use Prentice formulation
        nll_agg_obs(y,i) -= keep.cdf_lower(keep_ind(y,i)) * log(squeeze(pgengamma(obsvec(keep_ind(y,i)), a, b, Q, 0, 1)));
        nll_agg_obs(y,i) -= keep.cdf_upper(keep_ind(y,i)) * log(1.0 - squeeze(pgengamma(obsvec(keep_ind(y,i)), a, b, Q, 0, 1)));
      } else { //delta model
        Type p_pos = annual_nll_pars(i,y,3);
        if(obs_pos(y,i)) nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (log(p_pos) + dgengamma(obsvec(keep_ind(y,i)), a, b, Q, 1, 1));
        else nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * log(1-p_pos);
      }
    }
    if(agg_obs_model(i) == 6){ //Tweedie
      Type mu = annual_nll_pars(i,y,0);
      Type power = annual_nll_pars(i,y,1);
      Type phi = annual_nll_pars(i,y,2);
      nll_agg_obs(y,i) -= keep(keep_ind(y,i)) *  dtweedie(obsvec(keep_ind(y,i)), mu, phi, power, 1); 
    }
  }
  return nll_agg_obs;
}


template <class Type>
matrix<int> simulate_obs_pos(vector<int> agg_obs_model, matrix<int> obs_pos, array<Type> annual_nll_pars, matrix<int> use_obs){
  int n_y = obs_pos.rows();
  int n_series = obs_pos.cols();
  for(int i = 0; i < n_series; i++) for(int y = 0; y < n_y; y++) if(use_obs(y,i)) {
    if(agg_obs_model(i) > 6)  {
      Type p_pos = 1.0;
      if(agg_obs_model(i) < 11) p_pos = annual_nll_pars(i,y,2);
      else p_pos = annual_nll_pars(i,y,3);
      Type x = rbinom(Type(1.0),p_pos);
      if(x>0.5) obs_pos(y,i) = 1;
      else obs_pos(y,i) = 0;
    }
  }
  return obs_pos;
}


template <class Type>
matrix<Type> simulate_agg_obs(vector<int> agg_obs_model, matrix<Type> agg_obs, matrix<int> obs_pos,
  array<Type> annual_nll_pars, matrix<int> use_obs){
  int n_y = agg_obs.rows();
  int n_series = agg_obs.cols();
  matrix<Type> agg_obs_out = agg_obs;
  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++) if(use_obs(y,i)) {
    if(agg_obs_model(i) < 3){ //log-normal
      agg_obs(y,i) = exp(rnorm(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1)));
    }
    if((agg_obs_model(i) >2) & (agg_obs_model(i) < 5)){
      agg_obs(y,i) = rgamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1));
    }
    if(agg_obs_model(i) == 5){
      agg_obs(y,i) = rgengamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1), annual_nll_pars(i,y,2));
    }
    if(agg_obs_model(i) == 6){ //Tweedie
      agg_obs(y,i) = rTweedie(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1), annual_nll_pars(i,y,2));
    }
    if(agg_obs_model(i) > 6){
      if(obs_pos(y,i)) { //delta models, whether index positive has already been simulated
        if(agg_obs_model(i) < 9){ //delta log-normal
          agg_obs(y,i) = exp(rnorm(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1)));
        }
        if((agg_obs_model(i) > 8) & (agg_obs_model(i) < 11)){ //delta gamma
          agg_obs(y,i) = rgamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1));
        }
        if(agg_obs_model(i) == 11){ //delta gen gamma
          agg_obs(y,i) = rgengamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1),annual_nll_pars(i,y,2));
        }
      } else agg_obs(y,i) = 0.0;
    }
  }
  return agg_obs;
}

template <class Type>
vector<Type> sim_agg_obs_in_obsvec(vector<Type> obsvec, matrix<int> keep_ind, vector<int> agg_obs_model, matrix<Type> agg_obs, matrix<int> obs_pos, matrix<int> use_obs){
  int n_y = use_obs.rows();
  int n_series = use_obs.cols();
  vector<Type> obsvec_out = obsvec;
  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++){
    if(use_obs(y,i)) {
      obsvec_out(keep_ind(y,i)) = agg_obs(y,i);
      if(agg_obs_model(i) < 3) obsvec_out(keep_ind(y,i)) = log(agg_obs(y,i)); //normal for log-index
      if((agg_obs_model(i) == 7) | (agg_obs_model(i) == 8)) {
        if(obs_pos(y,i)) obsvec_out(keep_ind(y,i)) = log(agg_obs(y,i)); //normal for log of positive index 
        // else it is already set to 0
      }
    }
  }
  return obsvec_out;
}

template <class Type>
vector<Type> obsvec_to_paa(int i, int y, vector<Type> obsvec, vector<int> agesvec, matrix<int> use_paa, array<int> keep_paa, 
  vector<int> age_comp_model, int n_ages){
  vector<Type> paa(n_ages);
  paa.setZero();
  if(use_paa(y,i)) {
    vector<int> ages_obs_y = agesvec.segment(keep_paa(i,y,0), keep_paa(i,y,1));
    vector<Type> tf_paa_obs = obsvec.segment(keep_paa(i,y,0),keep_paa(i,y,1));
    paa = make_paa(tf_paa_obs, age_comp_model(i), ages_obs_y, paa);
  }
  return paa;
}
