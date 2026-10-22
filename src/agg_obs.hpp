
//find sigma given CV+ and nu for generalized gamma distribution (Rigby et al. 2019 parameterization)
template<class Type>
struct gengamma_trans_sigma_objective {
  /* Data and parameter objects for calculation: */
  Type nu;
  Type CV;
  int nu_pos=1;
  int trace=0;

  /* Constructor */
  gengamma_trans_sigma_objective(
  Type nu_,
  Type CV_,
  int nu_pos_,
  int trace_) :
    nu(nu_), CV(CV_), nu_pos(nu_pos_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> trans_sigma) {
    T nu_T = T(nu);
    T sigma = exp(trans_sigma(0));
    if(!nu_pos) sigma = 1.0/(2.0*pow(nu_T*nu_T,0.5)*(1.0+exp(-trans_sigma(0)))); //0<sigma < 1/(2*|nu|) for nu<0
    T theta = pow(sigma*nu_T, -2.0);
    vector<T> diff(1), diff2(1);
    diff(0) = T(CV)*T(CV) + 1.0 - exp(lgamma(theta) + lgamma(theta + 2.0/nu_T) - 2.0*lgamma(theta + 1.0/nu_T));
    diff2(0) = diff(0) * diff(0);
    if(trace) {
      see(T(CV)*T(CV) + 1.0);
      see(exp(lgamma(theta) + lgamma(theta + 2.0/T(nu)) - 2.0*lgamma(theta + 1.0/nu_T)));
      see(diff);
      see(theta);
      see(nu_T);
    }
    return diff2;
  }
};

//find nuplus given CV+ and sigma for generalized gamma distribution (Rigby et al. 2019 parameterization)
template<class Type>
struct gengamma_log_nuplus_objective {
  /* Data and parameter objects for calculation: */
  Type sigma;
  Type CV;
  int trace=0;

  /* Constructor */
  gengamma_log_nuplus_objective(
  Type sigma_,
  Type CV_,
  int trace_) :
    sigma(sigma_), CV(CV_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> log_nuplus) {
    T nuplus = exp(log_nuplus(0));
    T nu = nuplus - 1/(2*pow(T(sigma),2.0));
    T theta = pow(T(sigma)*nu, -2.0);
    vector<T> diff(1), diff2(1);
    diff(0) = T(CV)*T(CV) + 1.0 - exp(lgamma(theta) + lgamma(theta + 2.0/nu) - 2.0*lgamma(theta + 1.0/nu));
    diff2(0) = diff(0) * diff(0);
    if(trace) {
      see(T(CV)*T(CV) + 1.0);
      see(exp(lgamma(theta) +lgamma(theta+2.0/nu) - 2.0*lgamma(theta + 1.0/nu)));
      see(diff);
      see(theta);
      see(sigma);
    }
    return diff2;
  }
};

template <class Type>
matrix<Type> get_gengamma_trans_sigma(Type CV, Type nu, int n = 30, int nu_pos = 1, int trace = 0){
  //Rigby et al. 2019 parameterization: find sigma given CV and nu. When nu is negative, sigma has a narrow range of values, but it is just bounded below by 0 when nu is positive.
  matrix<Type> trans_sigma_iter(n,3);
  trans_sigma_iter.setZero();
  gengamma_trans_sigma_objective<Type> trans_sigma_obj(nu, CV, nu_pos, trace);
  vector<Type> trans_sigma_i(1);
  trans_sigma_i(0) = 0.0;
  for (int i=0; i<n-1; i++) {
    trans_sigma_i(0) = trans_sigma_iter(i,0); //current log_sigma
    trans_sigma_iter(i,1) = trans_sigma_obj(trans_sigma_i)(0); //f(log_sigma)
    trans_sigma_iter(i,2) = autodiff::jacobian(trans_sigma_obj,trans_sigma_i)(0,0); //f'(log_sigma)
    trans_sigma_iter(i+1,0) = trans_sigma_iter(i,0) - trans_sigma_iter(i,1)/trans_sigma_iter(i,2); //next log_sigma
  }
  return(trans_sigma_iter);
}

template <class Type>
matrix<Type> get_gengamma_log_nuplus(Type CV, Type sigma, int n = 30, int trace = 0){
  //Rigby et al. 2019 parameterization: find nu given CV and sigma. Generally, the relationship between CV and nu is not monotonic so 2 values of nu can give the same CV.
  //nu constrained to be greater than -1/(2*sigma^2) so that variance is defined. But pdf cannot be evaluated at 0.
  matrix<Type> log_nuplus_iter(n,3);
  log_nuplus_iter.setZero();
  gengamma_log_nuplus_objective<Type> log_nuplus_obj(sigma, CV, trace);
  vector<Type> log_nuplus_i(1);
  log_nuplus_i(0) = 0.0;
  for (int i=0; i<n-1; i++) {
    log_nuplus_i(0) = log_nuplus_iter(i,0); //current log_nuplus
    log_nuplus_iter(i,1) = log_nuplus_obj(log_nuplus_i)(0); //f(log_nuplus)
    log_nuplus_iter(i,2) = autodiff::jacobian(log_nuplus_obj,log_nuplus_i)(0,0); //f'(log_nuplus)
    log_nuplus_iter(i+1,0) = log_nuplus_iter(i,0) - log_nuplus_iter(i,1)/log_nuplus_iter(i,2); //next log_nuplus
  }
  return(log_nuplus_iter);
}

template <class Type>
Type dgengamma(Type x, Type mu, Type sigma, Type nu, int give_log=0){
  //Rigby et al. 2019 parameterization
  //Type nu = nu_in;
  //if(use_var_bound) nu -= 1/(2*sigma*sigma); //nu constrained to be greater than -1/(2*sigma^2) so that variance is defined. But pdf cannot be evaluated at 0.
  Type theta = pow(sigma*nu,-2.0); 
  Type z = pow(x/mu, nu);
  Type logres = 0.0;
  logres = 0.5*log(nu*nu) + theta*log(theta) + theta*log(z) - theta*z - lgamma(theta) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type pgengamma(Type x, Type mu, Type sigma, Type nu){
  //Rigby et al. 2019 parameterization
  Type res = 0.0;
  Type theta = pow(sigma*nu,-2.0); //= alpha shape parameter
  Type z = pow(x/mu, nu);
  res = pgamma(z,theta,1/theta); //TMB only uses shape and scale (not rate). E(z) = alpha * beta = 1, V(z) = (sigma*nu)^2
  return res;
}

template<class Type>
Type rgengamma(Type mu, Type sigma, Type nu){
  Type alpha = pow(sigma*nu,-2.0); 
  Type beta = 1/alpha;
  Type z = rgamma(alpha,beta);
  Type x = mu * pow(z, 1/nu);
  return  x;
}

template<class Type>
Type rgengamma_pred(Type pred, Type sigma, Type nu){
  Type alpha = pow(sigma*nu,-2.0); // alpha=theta
  Type beta = 1/alpha;
  Type mu = pred * exp(-lgamma(alpha + 1/nu) + lgamma(alpha) + log(alpha)/nu);
  Type z = rgamma(alpha,beta);
  Type x = mu * pow(z, 1/nu);
  return  x;
}

// Simulate from tweedie/CPG
// Adapted from tweedie::rtweedie function in R
template<class Type>
Type rTweedie( Type mu, Type phi, Type power){
 Type lambda = pow(mu, Type(2.0) - power) / (phi * (Type(2.0) - power));
 Type alpha = (Type(2.0) - power) / (power-Type(1.0)); //alpha>0 unlike tweedie::rtweedie
 Type gam = phi * (power - Type(1.0)) * pow(mu, power - Type(1.0));
 Type N = rpois(lambda);
 Type B = rgamma(N * alpha, gam);   /// Using Shape-Scale parameterization
 return B;
}

template<class Type>
struct CPG_log_lambda_objective_alpha {
  /* Data and parameter objects for calculation: */
  Type alpha;
  Type cvplus;
  int trace;

  /* Constructor */
  CPG_log_lambda_objective_alpha(
  Type alpha_,
  Type cvplus_,
  int trace_) :
    alpha(alpha_), cvplus(cvplus_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> log_lambda) { //find CPG lambda given CV+ and alpha.
    T lambda = exp(log_lambda(0));
    T p0 = exp(-lambda);
    T cv2hat = (1-p0)*(T(alpha) + 1)/(lambda*T(alpha))  - p0; //cannot be < 0
    vector<T> diff(1);
    diff(0) = cv2hat - T(cvplus)*T(cvplus);
    if(trace) see(cv2hat);
    return diff * diff;
  }
};

template<class Type>
struct CPG_log_lambda_objective_beta {
  /* Data and parameter objects for calculation: */
  Type beta;
  Type cvplus;
  Type muplus;
  Type lambda_min;
  int trace;

  /* Constructor */
  CPG_log_lambda_objective_beta(
  Type beta_,
  Type cvplus_,
  Type muplus_,
  Type lambda_min_,
  int trace_) :
    beta(beta_), cvplus(cvplus_), muplus(muplus_), lambda_min(lambda_min_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> log_lambda) { //find CPG lambda given CV+, E+, beta, and minimum lambda allowed by CV+.
    T lambda = exp(log_lambda(0))+T(lambda_min);
    T p0 = exp(-lambda);
    T cv2hat = (1-p0)/lambda + T(beta)/T(muplus) - p0; //cannot be < 0
    if(trace) see(lambda);
    if(trace) see(p0);
    if(trace) see(beta);
    if(trace) see(muplus);
    if(trace) see(cvplus);
    if(trace) see(cv2hat);
    vector<T> diff(1);
    diff(0) = cv2hat - T(cvplus)*T(cvplus);
    if(trace) see(diff);
    if(trace) see(diff*diff);
    return diff * diff;
  }
};

template <class Type>
matrix<Type> get_CPG_log_lambda(int given_alpha, Type alpha, Type beta, Type lambda_min, Type muplus, Type cvplus, int n = 30, int trace = 0){
  matrix<Type> log_lambda_iter(n,3);
  log_lambda_iter.setZero();
  vector<Type> log_lambda_i(1);
  log_lambda_iter(0,0) = log(4.0); //start on the right side of the peak
  if(given_alpha) { //lambda = f(CV+, alpha)
    CPG_log_lambda_objective_alpha<Type> lambda_obj(alpha, cvplus, trace);
    for (int i=0; i<n-1; i++) {
      log_lambda_i(0) = log_lambda_iter(i,0);
      log_lambda_iter(i,1) = lambda_obj(log_lambda_i)(0); //f(log_lambda)
      log_lambda_iter(i,2) = autodiff::jacobian(lambda_obj,log_lambda_i)(0,0); //f'(log_lambda)
      log_lambda_iter(i+1,0) = log_lambda_iter(i,0) - log_lambda_iter(i,1)/log_lambda_iter(i,2); 
    }
  }
  else { //lambda = f(CV+, E+, beta, lambda_min);
    CPG_log_lambda_objective_beta<Type> lambda_obj(beta, cvplus, muplus, lambda_min, trace);
    for (int i=0; i<n-1; i++) {
      log_lambda_i(0) = log_lambda_iter(i,0);
      log_lambda_iter(i,1) = lambda_obj(log_lambda_i)(0); //f(log_lambda)
      log_lambda_iter(i,2) = autodiff::jacobian(lambda_obj,log_lambda_i)(0,0); //f'(log_lambda)
      log_lambda_iter(i+1,0) = log_lambda_iter(i,0) - log_lambda_iter(i,1)/log_lambda_iter(i,2); 
    }
  }
  return(log_lambda_iter);
}

template<class Type>
struct CPG_lambda_min_objective {
  /* Data and parameter objects for calculation: */
  Type cvplus;
  Type lambda_at_max;
  int trace;

  /* Constructor */
  CPG_lambda_min_objective(
  Type cvplus_,
  Type lambda_at_max_,
  int trace_) :
    cvplus(cvplus_), lambda_at_max(lambda_at_max_), trace(trace_) {}

  template <typename T> //This allows you to differentiate the function wrt whatever is after operator() on line below
  vector<T> operator()(vector<T> lambda) { //find CPG minimum lambda consistent with given CV+ when CV+<0.55.
    T lambdaT = T(lambda_at_max) + lambda(0); //use the right side of the peak for this non-monotonic function
    T p0 = exp(-lambdaT);
    T cv2minhat = (1-p0)/lambdaT - p0; //cannot be < 0
    vector<T> diff(1);
    diff(0) = T(cvplus)*T(cvplus) - cv2minhat;
    if(trace) see(cv2minhat);
    return diff * diff;
  }
};

template <class Type>
matrix<Type> get_CPG_lambda_min(Type cvplus, int n = 30, int trace = 0){
  matrix<Type> lambda_min_iter(n,3);
  lambda_min_iter.setZero();
  CPG_lambda_min_objective<Type> lambda_min_obj(cvplus, 1.79328, trace);
  vector<Type> lambda_min_i(1); //this is lambda_min - 1.79328. the function adds 1.79328
  lambda_min_iter(0,0) = 2.0; //seems to be a good starting point
  for (int i=0; i<n-1; i++) {
    lambda_min_i(0) = lambda_min_iter(i,0);
    lambda_min_iter(i,1) = lambda_min_obj(lambda_min_i)(0); //f(lambda_min)
    lambda_min_iter(i,2) = autodiff::jacobian(lambda_min_obj,lambda_min_i)(0,0); //f'(lambda_min)
    lambda_min_iter(i+1,0) = lambda_min_iter(i,0) - lambda_min_iter(i,1)/lambda_min_iter(i,2); 
  }
  return(lambda_min_iter); //lambda-1.79328
}

template <class Type>
matrix<Type> get_lambda_min(matrix<Type> cvplus, vector<int> catch_model, int n_iter = 30, int trace = 0){
  int n_series = cvplus.cols();
  int n_years = cvplus.rows();
  matrix<Type> lambda_min(n_years, n_series);
  lambda_min.setZero();
  for(int i = 0; i < n_series; i++) {
    if(catch_model(i) == 7) for(int y = 0; y < n_years; y++) lambda_min(y,i) = 1.79328;
    if(catch_model(i) == 8) for(int y = 0; y < n_years; y++){
        matrix<Type> lambda_min_iter = get_CPG_lambda_min(cvplus(y,i), n_iter, trace);
        lambda_min(y,i) = lambda_min_iter(n_iter-1,0);
        lambda_min(y,i) = 1.79328 + lambda_min(y,i)/(1 + exp(1000*(cvplus(y,i)-0.5463))); //if cv>0.55 minimum lambda is 0.
        lambda_min(y,i) = lambda_min(y,i)/(1 + exp(1000*(cvplus(y,i)-0.5463))); //if cv>0.55 minimum lambda is 0.
    }
  }
  return lambda_min;
}

template <class Type>
array<Type> get_annual_agg_obs_pars(vector<int> agg_obs_model, matrix<Type> pred, matrix<Type> obs_cv, matrix<Type> effort, matrix<Type> agg_obs_pars, matrix<Type> lambda_min, int bias_correct_oe, int n_iter = 30, int trace = 0){
/*
     1: lognormal-1
     2: lognormal-2
     3: gamma-1
     4: gamma-2
     5: gengamma
     6: CPG-1
     7: CPG-2
     8: CPG-3
     9: delta-lognormal-1
    10: delta-lognormal-2
    11: delta-gamma-1
    12: delta-gamma-2
    13: delta-gengamma
*/

  int n_y = pred.rows();
  int n_series = pred.cols();
  array<Type> annual_nll_pars(n_series,n_y,4);
  annual_nll_pars.setZero();

  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++) {
    if((agg_obs_model(i) == 1) | (agg_obs_model(i) == 9)){ //normal for log-obs: up to 1 free parameters (scalar)  (delta: up to 3)
      annual_nll_pars(i,y,0) = log(pred(y,i)); //mu for normal pdf of log-obs
      Type obs_sd = pow(log(pow(obs_cv(y,i), 2.0) + 1), 0.5) * exp(agg_obs_pars(i,0)); //SD = sigma * phi
      if(bias_correct_oe) annual_nll_pars(i,y,0) -= 0.5*pow(obs_sd,2.0);
      annual_nll_pars(i,y,1) = obs_sd; //sigma for normal pdf of log-obs
    }
    if((agg_obs_model(i) == 2) | (agg_obs_model(i) == 10)){ //normal for log-obs: with CV and SD(log(catch)) required to be greater than input  (delta: up to 3)
      annual_nll_pars(i,y,0) = log(pred(y,i)); //mu for normal pdf of log-obs
      Type obs_sd = pow(log(pow(obs_cv(y,i), 2.0) + 1), 0.5) * (1.0 + exp(agg_obs_pars(i,0))); //SD = sigma * (1+phi)
      if(bias_correct_oe) annual_nll_pars(i,y,0) -= 0.5*pow(obs_sd,2.0);
      annual_nll_pars(i,y,1) = obs_sd; //sigma for normal pdf of log-obs
    }
    if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 11)){ //gamma: up to 1 free parameters (scalar) (delta: up to 3)
      //gamma shape and scale parameters are returned. CV is input CV multiplied by scalar CV = cv_in * scalar
      // var(x) = E(x)^2 * CV^2 = (alpha * beta)^2 * CV^2
      //Type cv = exp(agg_obs_pars(i,0)) * obs_sd(y,i)/obs(y,i);
      Type cv = obs_cv(y,i) * exp(agg_obs_pars(i,0));
      annual_nll_pars(i,y,0) = 1/(cv * cv); //gamma alpha (shape) parameter:  mu/beta
      annual_nll_pars(i,y,1) = pred(y,i) * cv * cv; //gamma beta(scale) parameter beta = mu/alpha
    }
    if((agg_obs_model(i) == 4) | (agg_obs_model(i) == 12)){ //gamma:  with variance required to be > alpha * beta^2  (delta: up to 3)
      //gamma shape and scale parameters are returned. CV must be greater than input: CV = cv_in * (scalar + 1)
      // var(x) = E(x)^2 * CV^2 = (alpha * beta)^2 * CV^2
      //Type cv = (exp(agg_obs_pars(i,0)) + 1.0) * obs_sd(y,i)/obs(y,i);
      Type cv = obs_cv(y,i) * (exp(agg_obs_pars(i,0)) + 1.0);
      annual_nll_pars(i,y,0) = 1/(cv * cv); //gamma alpha (shape) parameter:  mu/beta
      annual_nll_pars(i,y,1) = pred(y,i) * cv * cv; //gamma beta(scale) parameter beta = mu/alpha
    }
    if((agg_obs_model(i) == 5) | (agg_obs_model(i) == 13)){ //generalized gamma: default 1 free parameters (delta: up to 3)
      //Type cv = pow(CV2,0.5);
      //Type cv = obs_sd(y,i) /obs(y,i);
      Type nu = exp(agg_obs_pars(i,0));
      matrix<Type> log_sigma_iter = get_gengamma_trans_sigma(obs_cv(y,i),nu,n_iter,1,trace); //iterations for finding log(sigma) of Generalized Gamma RVs. Last argument assumes nu is strictly positive, so that sigma > 0
      Type sigma = exp(log_sigma_iter(n_iter-1,0));
      Type theta = pow(nu*sigma,-2.0);
      Type mu = pred(y,i) * exp(-lgamma(theta + 1/nu) + lgamma(theta) + log(theta)/nu);
      annual_nll_pars(i,y,0) = mu; //scale par
      annual_nll_pars(i,y,1) = sigma; //shape par
      annual_nll_pars(i,y,2) = nu;
    }
    if(agg_obs_model(i) > 8){ //delta-models: up to 3 free parameter
      //9: delta-lognormal, 10: delta-lognormal(scalar>1), 11: delta-gamma, 12: delta-gamma(var>mu/beta), 13: delta-generalized gamma
      Type beta0 = agg_obs_pars(i,1); //intercept/constant lambda parameter
      Type beta1 = exp(agg_obs_pars(i,2)); //rate of positive obs with pred_obs (q*sum(sel*N)))
      Type beta2 = exp(agg_obs_pars(i,3)); //rate of positive obs with measure of "effort")
      Type lambda = exp(beta0 + beta1*log(pred(y,i)) + beta2 * log(effort(y,i)));
      Type p_pos = (1 - exp(-lambda)); //probability of positive obs: 1 - prob Poisson RV = 0 for Tweedie.
      annual_nll_pars(i,y,3) = p_pos; //probablity of positive obs
    }
    if((agg_obs_model(i) > 5) & (agg_obs_model(i) < 9)){ //Tweedie (compound Poisson-Gamma): 1-3 free parameters
      Type cv_pos = obs_cv(y,i); //obs_sd(y,i) /obs(y,i);
      Type mu_pos = pred(y,i);
      Type alpha = 0, lambda = 0, beta = 0, p_pos = 0;
      
      //for model == 6, alpha free, lambda = f(alpha,CV), beta = f(alpha,lambda,pred), but 0<alpha<1 to ensure unique lambda
      if(agg_obs_model(i) == 6){ 
        alpha = 1/(1+ exp(-agg_obs_pars(i,0))); // 0<alpha<1 to ensure one value of lambda can be found in get_CPG_log_lambda()
        //iterations for finding log(lambda) of CPG (mean of Poisson for number of iid gammma RVs).
        matrix<Type> log_lambda_iter = get_CPG_log_lambda(1, alpha, Type(0), Type(0), mu_pos, cv_pos, n_iter, trace); 
        lambda = exp(log_lambda_iter(n_iter-1,0)); //cannot be < 0, lower limit is not required
        p_pos = 1.0 - exp(-lambda);
        beta = p_pos * mu_pos/(lambda * alpha); //cannot be < 0
      }

      //for model == 7, beta constant fraction of pred*CV^2 (free, but 0 < beta < pred*CV^2), lambda = f(CV, beta, lambda_min), solved, alpha = f(lambda, pred, beta)
      if(agg_obs_model(i) == 7){ 
        //for a given CV, there is a minimum lambda so that alpha and beta > 0
        Type beta_max = mu_pos * cv_pos * cv_pos;
        beta = beta_max/(1 + exp(-agg_obs_pars(i,0)));
        //iterations for finding log(lambda) of CPG (mean of Poisson for number of iid gammma RVs).
        matrix<Type> log_lambda_iter = get_CPG_log_lambda(0, Type(0), beta, lambda_min(y,i), mu_pos, cv_pos, n_iter, trace);
        lambda = lambda_min(y,i) + exp(log_lambda_iter(n_iter-1,0)); //cannot be < 0, lower limit is not required
        p_pos = 1.0 - exp(-lambda);
        alpha = p_pos * mu_pos/(lambda*beta);
      }

      //for model == 8, lambda freely estimated possibly a function of effort and mu_pos, alpha = f(lambda, CV), beta = f(alpha,lambda,pred), but lower bound of lambda defined by CV to ensure alpha>0
      if(agg_obs_model(i) == 8){ 
        //for a given CV, there is a minimum lambda so that alpha and beta > 0
        Type lam0 = agg_obs_pars(i,0); //constant or intercept log(lambda) parameter 
        Type lam1 = exp(agg_obs_pars(i,1)); //>0 power for pred_obs (q*sum(sel*N))) effect on lambda, increased lambda = decreased p_0
        Type lam2 = exp(agg_obs_pars(i,2)); //>0 power for "effort" (q*sum(sel*N))) effect on lambda.
        lambda = lambda_min(y,i) + exp(lam0 + lam1*log(pred(y,i)) + lam2 * log(effort(y,i)));
        p_pos = (1 - exp(-lambda)); //probability of positive obs: 1 - prob Poisson RV = 0 for CPG.
        Type part = (cv_pos*cv_pos - p_pos/lambda + 1-p_pos);
        beta = mu_pos * part;
        alpha = p_pos/(lambda*part);
      }

      Type mu = p_pos * mu_pos; //mean of Tweedie/CPG
      Type power = (2.0 + alpha)/(1.0 + alpha); //Tweedie parameter: 1<power<2 for CPG
      Type phi = (alpha + 1.0)*beta*pow(mu, -1.0/(alpha+1)); //Tweedie parameter: V(x) = phi * mu^power
      
      annual_nll_pars(i,y,0) = mu;
      annual_nll_pars(i,y,1) = phi;
      annual_nll_pars(i,y,2) = power;
    }
  }
  return annual_nll_pars;
}


template <class Type>
matrix<Type> get_nll_agg_obs(vector<int> agg_obs_model, matrix<int> obs_pos, array<Type> annual_nll_pars,
  vector<Type> obsvec, matrix<int> use_obs, matrix<int> keep_ind, data_indicator<vector<Type>, Type> keep){
  int n_series = annual_nll_pars.dim(0);
  int n_y = annual_nll_pars.dim(1);

/*
   1: lognormal-1
     2: lognormal-2
     3: gamma-1
     4: gamma-2
     5: gengamma
     6: CPG-1
     7: CPG-2
     8: CPG-3
     9: delta-lognormal-1
    10: delta-lognormal-2
    11: delta-gamma-1
    12: delta-gamma-2
    13: delta-gengamma
*/

  //for delta(hurdle) models and Tweedie, need to make sure call to TMB::oneStepPredict has discreteSupport = 0 and range = c(0,Inf)
  matrix<Type> nll_agg_obs(n_y,n_series);
  nll_agg_obs.setZero();

  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++) if(use_obs(y,i)){
    if((agg_obs_model(i) <3) | (agg_obs_model(i) == 9) | (agg_obs_model(i) == 10)) { //(delta) normal for log-obs, two types
      Type mu = annual_nll_pars(i,y,0);
      Type sig = annual_nll_pars(i,y,1);
      if(agg_obs_model(i) < 3){
        nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (dnorm(obsvec(keep_ind(y,i)), mu, sig,1) - obsvec(keep_ind(y,i))); //log-normal
        nll_agg_obs(y,i) -= keep.cdf_lower(keep_ind(y,i)) * log(squeeze(pnorm(obsvec(keep_ind(y,i)), mu, sig)));
        nll_agg_obs(y,i) -= keep.cdf_upper(keep_ind(y,i)) * log(1.0 - squeeze(pnorm(obsvec(keep_ind(y,i)), mu, sig)));
      } else { //delta model
        Type p_pos = annual_nll_pars(i,y,3);
        if(obs_pos(y,i)) nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (log(p_pos) + dnorm(obsvec(keep_ind(y,i)), mu, sig,1) - obsvec(keep_ind(y,i))); //log-normal
        else nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * log(1-p_pos);
      }
    }
    if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 4) | (agg_obs_model(i) == 11) | (agg_obs_model(i) == 12)) {//(delta) gamma, two types
      Type alpha = annual_nll_pars(i,y,0);
      Type beta = annual_nll_pars(i,y,1); //scale, not rate
      if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 4)){ //not delta model
        nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * dgamma(obsvec(keep_ind(y,i)), alpha, beta,1);
        nll_agg_obs(y,i) -= keep.cdf_lower(keep_ind(y,i)) * log(squeeze(pgamma(obsvec(keep_ind(y,i)), alpha, beta)));
        nll_agg_obs(y,i) -= keep.cdf_upper(keep_ind(y,i)) * log(1.0 - squeeze(pgamma(obsvec(keep_ind(y,i)), alpha, beta)));
      } else { //delta model
        Type p_pos = annual_nll_pars(i,y,3);
        if(obs_pos(y,i)) nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (log(p_pos) + dgamma(obsvec(keep_ind(y,i)), alpha, beta,1));
        else nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * log(1-p_pos);
      }
    }
    if((agg_obs_model(i) == 5) | (agg_obs_model(i) == 13)){ //(delta) generalized gamma
      Type mu = annual_nll_pars(i,y,0); //mean par
      Type sigma = annual_nll_pars(i,y,1); //shape par
      Type nu = annual_nll_pars(i,y,2); //could be configured to be negative, but complicated.
      if(agg_obs_model(i) == 5){
        nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * dgengamma(obsvec(keep_ind(y,i)), mu, sigma, nu, 1); //Use Rigby et al. (2019) formulation
        nll_agg_obs(y,i) -= keep.cdf_lower(keep_ind(y,i)) * log(squeeze(pgengamma(obsvec(keep_ind(y,i)), mu, sigma, nu)));
        nll_agg_obs(y,i) -= keep.cdf_upper(keep_ind(y,i)) * log(1.0 - squeeze(pgengamma(obsvec(keep_ind(y,i)), mu, sigma, nu)));
      } else { //delta model
        Type p_pos = annual_nll_pars(i,y,3);
        if(obs_pos(y,i)) nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * (log(p_pos) + dgengamma(obsvec(keep_ind(y,i)), mu, sigma, nu, 1));
        else nll_agg_obs(y,i) -= keep(keep_ind(y,i)) * log(1-p_pos);
      }
    }
    if((agg_obs_model(i) > 5) & (agg_obs_model(i) < 9)){ //Tweedie ll, 12, 13 with different CPG parameters held constant
      Type mu = annual_nll_pars(i,y,0);
      Type phi = annual_nll_pars(i,y,1);
      Type power = annual_nll_pars(i,y,2);
      nll_agg_obs(y,i) -= keep(keep_ind(y,i)) *  dtweedie(obsvec(keep_ind(y,i)), mu, phi, power, 1); 
    }
  }
  return nll_agg_obs;
}


//just for delta models (not CPG)
template <class Type>
matrix<int> simulate_obs_pos(vector<int> agg_obs_model, matrix<int> obs_pos, array<Type> annual_nll_pars, matrix<int> use_obs){
  int n_y = obs_pos.rows();
  int n_series = obs_pos.cols();
  for(int i = 0; i < n_series; i++) if(agg_obs_model(i) > 8) for(int y = 0; y < n_y; y++) if(use_obs(y,i)) {
    obs_pos(y,i) = CppAD::Integer(rbinom(Type(1.0),annual_nll_pars(i,y,3)));
  }
  return obs_pos;
}


template <class Type>
matrix<Type> simulate_agg_obs(vector<int> agg_obs_model, matrix<Type> agg_obs, matrix<int> obs_pos,
  array<Type> annual_nll_pars, matrix<int> use_obs){
  int n_y = agg_obs.rows();
  int n_series = agg_obs.cols();
  for(int i = 0; i < n_series; i++) for(int y = 0; y < n_y; y++) if(use_obs(y,i)) {
    if(agg_obs_model(i) < 3){ //log-normal
      agg_obs(y,i) = exp(rnorm(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1)));
    }
    if((agg_obs_model(i) == 3) | (agg_obs_model(i) == 4)){
      agg_obs(y,i) = rgamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1));
    }
    if(agg_obs_model(i) == 5){
      agg_obs(y,i) = rgengamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1), annual_nll_pars(i,y,2));
    }
    if(agg_obs_model(i) > 8){
      if(obs_pos(y,i)) { //delta models, whether index is positive has already been simulated
        if((agg_obs_model(i) == 9) | (agg_obs_model(i) == 10)){ //delta log-normal
          agg_obs(y,i) = exp(rnorm(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1)));
        }
        if((agg_obs_model(i) == 11) | (agg_obs_model(i) == 12)){ //delta gamma
          agg_obs(y,i) = rgamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1));
        }
        if(agg_obs_model(i) == 13){ //delta gen gamma
          agg_obs(y,i) = rgengamma(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1),annual_nll_pars(i,y,2));
        }
      } else agg_obs(y,i) = 0.0;
    }
    if((agg_obs_model(i) > 5) & (agg_obs_model(i) < 9)){ //Tweedie
      agg_obs(y,i) = rTweedie(annual_nll_pars(i,y,0), annual_nll_pars(i,y,1), annual_nll_pars(i,y,2));
    }
  }
  return agg_obs;
}

template <class Type>
matrix<int> set_obs_pos_tweedie(vector<int> agg_obs_model, matrix<int> obs_pos, matrix<Type> agg_obs, matrix<int> use_obs) {
  int n_y = use_obs.rows();
  int n_series = use_obs.cols();
  for(int i = 0; i < n_series; i++) if((agg_obs_model(i) > 5) & (agg_obs_model(i) < 9)){
    for(int y = 0; y < n_y; y++) if(use_obs(y,i)){
      obs_pos(y,i) = 1;
      if(agg_obs(y,i)<1e-10) obs_pos(y,i) = 0;
    }
  }
  return obs_pos;
}

template <class Type>
vector<Type> sim_agg_obs_in_obsvec(vector<Type> obsvec, matrix<int> keep_ind, matrix<Type> agg_obs, matrix<int> obs_pos, matrix<int> use_obs, vector<int> agg_obs_model){
  int n_y = use_obs.rows();
  int n_series = use_obs.cols();
  vector<Type> obsvec_out = obsvec;
  for(int y = 0; y < n_y; y++) for(int i = 0; i < n_series; i++){
    if(use_obs(y,i)) {
      if((agg_obs_model(i) < 3) | (((agg_obs_model(i) == 9) | (agg_obs_model(i) == 10)) & obs_pos(y,i))) obsvec_out(keep_ind(y,i)) = log(agg_obs(y,i)); //lognormal
      else obsvec_out(keep_ind(y,i)) = agg_obs(y,i);
    }
  }
  return obsvec_out;
}

template <class Type>
matrix<Type> get_cvs(matrix<Type> cv_in, vector<Type> cv_re, matrix<int> obs_pos, matrix<int> use_obs, matrix<Type> cv_repars) {

  int n_series = cv_in.cols();
  int n_y = cv_in.rows();
  matrix<Type> cvs = cv_in;
  int cnt = 0;
  int n_re = cv_re.size();
  if(n_re>0) for(int i = 0; i < n_series; i++) for(int y = 0; y < n_y; y++) {
    if(use_obs(y,i)) {
      if(obs_pos(y,i)){
        cvs(y,i) = cv_in(y,i);
      } else{
        cvs(y,i) = exp(cv_re(cnt));
        cnt++;
      }
    }
  }
  return(cvs);
}

template <class Type>
vector<Type> simulate_cv_re(vector<int> cv_use_re, matrix<Type> cvs, vector<Type> cv_re, matrix<Type> cv_repars, matrix<int> use_obs, matrix<int> obs_pos) {
  int n_series = cv_use_re.size();
  int n_y = cvs.rows();
  vector<Type> cv_re_out = cv_re;
  cv_re_out.setZero();
  int cnt = 0;
  for(int i = 0; i < n_series;i++) if(cv_use_re(i)){
    Type mu = cv_repars(i,0);
    Type sigma = exp(cv_repars(i,1));
    Type rho = -1 + 2/(1 + exp(-cv_repars(i,2)));
    for(int y = 0; y<n_y; y++) if(use_obs(y,i) & (obs_pos(y,i)==0)){
      if(y == 0) {
        cv_re_out(cnt) = rnorm(mu, sigma*pow(1- pow(rho,2), -0.5));
      } else{
        cv_re_out(cnt) = rnorm(mu*(1-rho) + rho*log(cvs(y-1,i)),sigma);
      }
      cvs(y,i) = exp(cv_re_out(cnt));
      cnt++;
    }
  }
  return(cv_re_out);
}

template <class Type>
vector<Type> get_nll_cv_re(vector<int> cv_use_re, matrix<Type> cvs, matrix<Type> cv_repars, matrix<int> obs_pos){

  int n_series = cv_use_re.size();
  int n_y = cvs.rows();
  vector<Type>nll(n_series);
  nll.setZero();
  for(int i = 0; i < n_series;i++) if(cv_use_re(i)){
    Type mu = cv_repars(i,0);
    Type sigma = exp(cv_repars(i,1));
    Type rho = -1 + 2/(1 + exp(-cv_repars(i,2)));
    Type marg_sd = sigma*pow(1- pow(rho,2), -0.5);
    nll(i) -= dnorm(log(cvs(0,i)), mu,marg_sd,1);
    for(int y = 1; y<n_y; y++) nll(i) -= dnorm(log(cvs(y,i)), mu*(1-rho) + rho*log(cvs(y-1,i)), sigma, 1);
  }
  return(nll);
}

template <class Type>
matrix<Type> get_log_pred(matrix<Type> pred, array<Type> annual_nll_pars, vector<int> obs_model){
  int n_y = pred.rows();
  int n_series = pred.cols();
  matrix<Type> log_pred(n_y,n_series);
  log_pred.setZero();
  for(int i = 0; i < n_series; i++) for(int y = 0; y < n_y; y++) {
    if((obs_model(i) < 3) | (obs_model(i) == 9) | (obs_model(i) == 10)) {
      log_pred(y,i) = annual_nll_pars(i,y,0);
    } else {
      log_pred(y,i) = log(pred(y,i));
    }
  }
  return log_pred;
}

