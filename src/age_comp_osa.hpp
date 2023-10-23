//from Trijoulet et al. (2023). OSA residuals for composition data
//osa residual versions of dmultinom, ddirmultinom, dlogisticnormal
//see repo at 

template<class Type>
Type pbinom(Type x, Type n, Type p) {
  return 1. - pbeta(p, x+1, n-x);
}

//function written by Anders Nielsen to reorder the categories which is needed for OSA residuals of vector-valued observations
// that do not have likelihoods defined in TMB.
template<class Type>
vector<int> order(vector<Type> k){
  int n=k.size();
  vector<int> o(n);
  o.setZero();
  int at=-1;
  for(int i=0; i<n;++i){
    if(k(i)>0.5){o(++at) = i;}  
  }
  at=n;  
  for(int i=n-1; i>=0;--i){
    if(k(i)<0.5){o(--at) = i;}  
  }
  return o;
}

template <class Type>
Type dmultinom(vector<Type> x, vector<Type> p, vector<int> ages, data_indicator<vector<Type>, Type> keep, int give_log, int do_osa)
{
  Type logres = 0;
  vector<Type> p_x(ages.size());
  for(int i = 0; i < ages.size(); i++) p_x(i) = p(ages(i)-1);
  p_x /= sum(p_x);

  if(do_osa){
    vector<Type> k=keep;
    //vector<Type> l=keep.cdf_lower;
    //vector<Type> h=keep.cdf_upper;
    vector<int> o=order(k);
    x=x(o); p_x=p_x(o); k=k(o); //l=l(o); h=h(o);
    Type nUnused=asDouble(sum(x));
    Type pUsed=0;
    //Type cdf;
    for(int i=0; i<x.size(); ++i){
      if(i!=(x.size()-1)){
        vector<Type> x2(2), p2(2);
        //Type p_i = squeeze(p_x(i));
        //Type one_minus_pUsed_i = squeeze(1.0-pUsed);
        x2(0) = x(i);
        x2(1) = nUnused-x(i);
        p2(0) = squeeze(p_x(i));
        p2(0) /= squeeze((Type(1)-pUsed));
        p2(0) = squeeze(p2(0));
        // p2(0) = squeeze(p(i)/(Type(1)-pUsed));//(Type(1)-pUsed_i); //for log of any p = 0
        p2(1) = 1. - p2(0);
        logres += k(i) * dmultinom(x2,p2,1); //binomial the hard way.
        //logres += k(i)*dbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true);
        //cdf = pbinom(x(i),nUnused,p(i)/(Type(1)-pUsed));
        nUnused -= x(i);
        pUsed += p_x(i);
      }else{ // last index 
        logres += k(i)*Type(0);
        //cdf = Type(1);
      }
      //cdf = squeeze(cdf);
      //logres += l[i] * log( cdf );       // NaN protected
      //logres += h[i] * log( 1.0 - cdf ); // NaN protected
    }
  } else {
    logres = dmultinom(x,p_x,1);
  }
  if(give_log){
    return logres;
  }else{ 
    return exp(logres);
  }
}


template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> alpha, int do_log)
{
  Type phi = alpha.sum();
  int n = x.size();
  Type ll = lgamma(phi);
  for(int i = 0; i < n; i++) ll +=  -lgamma(alpha(i)) + (alpha(i) - 1.0) * log(x(i));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> p, Type phi, vector<int> ages, data_indicator<vector<Type>, Type> keep, int pool0, int give_log, int do_osa)
{
  vector<Type> p_pos(ages.size());
  p_pos.setZero();
  if(pool0){ //pooling zeros
    int index = 0;
    for(int a = 0; a < p.size(); a++)
    {
      p_pos(index) += p(a);
      if((a == ages(index)-1) & (index < (ages.size() - 1))) index++;
      //if(x(a) > Type(1.0e-15) & index < npos-1) index++;
    }
  } else { //missing zeros
    for(int a = 0; a < ages.size(); a++)
    {
      p_pos(a) = p(ages(a)-1);
      //x_pos = x(pos_ind);
    }
  }
  p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.

  vector<Type> alpha = p_pos * phi;
  Type logres = 0;
  if(do_osa){
    vector<Type> k=keep;
    //vector<Type> l=keep.cdf_lower;
    //vector<Type> h=keep.cdf_upper;
    vector<int> o=order(k);
    x=x(o); alpha=alpha(o); k=k(o); //l=l(o); h=h(o);
    
    int n = alpha.size();
    //Type cdf;
    Type sx = 1; // was: x.sum();
    Type sa = alpha.sum();
    sa -= alpha(0);
    logres = k(0)*dbeta(x(0),alpha(0),sa,true);
    //cdf = pbeta(x(0),alpha(0),sa);
    //cdf = squeeze(cdf);
    //logres += l(0) * log( cdf );       
    //logres += h(0) * log( 1.0 - cdf ); 
    
    for(int i=1; i<(n-1); ++i){
      sx -= x(i-1);
      sa -= alpha(i);
      logres += k(i)*(dbeta(x(i)/sx,alpha(i),sa,true)-log(sx));
      //cdf = pbeta(x(i)/sx,alpha(i),sa);
      //cdf = squeeze(cdf);
      //logres += l(i) * log( cdf );       
      //logres += h(i) * log( 1.0 - cdf ); 
    }
    logres += k(n-1)*Type(0);
    //cdf=Type(1);
    //cdf = squeeze(cdf);
    //logres += l(n-1) * log( cdf );       
    //logres += h(n-1) * log( 1.0 - cdf ); 
  } else {
    logres = ddirichlet(x,alpha,1);
  }
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}

//the usual D-M
template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type phi=sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
Type pbetabinom(Type x, Type N, Type alpha, Type beta, int do_log)
{
  //just sum up the probabilities at i = 0,...,x. Nothing fancy.
  int x_int = CppAD::Integer(x);
  vector<Type> p(2), aa(2);
  aa[0] = alpha;
  aa[1] = beta;
  Type Fx = 0.0;
  vector<Type> obs(2);
  for(int i = 0; i <= x_int; i++){
    obs[0] = i;
    obs[1] = N-i;
    Fx += ddirmultinom(obs, aa, 0); //dbetabinomial, just two categories.
  }
  if(do_log == 1) return log(Fx);
  else return Fx;
}

//the D-M as a series of conditional beta-binomials and added args for osa residuals
template<class Type> 
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, vector<int> ages, data_indicator<vector<Type>, Type> keep, int do_log, int do_osa)
{
  vector<Type> alpha_obs(ages.size());
  for(int i = 0; i < ages.size(); i++) alpha_obs(i) = alpha(ages(i)-1);
  
  Type ll = 0.0;
  if(do_osa){
    vector<Type> k=keep;
    //vector<Type> l=keep.cdf_lower;
    //vector<Type> h=keep.cdf_upper;
    
    vector<int> o=order(k);
    obs=obs(o); alpha_obs=alpha_obs(o); k=k(o); //l=l(o); h=h(o);   
    int dim = obs.size();
    vector<Type> alphas_a(2), obs_a(2);
    Type alp_sum = alpha_obs.sum();
    Type obs_sum = asDouble(obs.sum());
    for(int a = 0; a < dim-1; a++){
      obs_sum -= obs[a];
      alp_sum -= alpha_obs[a];
      obs_a(0) = obs(a);
      obs_a(1) = obs_sum;
      alphas_a(0) = alpha_obs(a);
      alphas_a(1) = alp_sum;
      ll += k(a) * ddirmultinom(obs_a, alphas_a, 1); //beta-binomial, just two categories
      //Type cdf = pbetabinom(obs_a(0), obs_a.sum(), alphas_a(0), alphas_a(1), 0);
      //cdf = squeeze(cdf);
      //ll += l(a) * log( cdf );
      //ll += h(a) * log( 1.0 - cdf );
    }
  }
  else{
    ll = ddirmultinom(obs,alpha_obs,1);
  }
  if(do_log) return ll;
  else return exp(ll);
}

template<class Type>
Type dmvnorm(vector<Type> x, vector<Type> mu,  matrix<Type> S, int do_log)
{
  using namespace density;
  MVNORM_t<Type> mvnorm(S);
  Type nll = mvnorm(x-mu);
  if(do_log) return -nll;
  else return exp(-nll);
}


//function to generate covariance matrix for logistic normal
template<class Type>
matrix<Type> make_AR_Sigma(int dim, vector<Type> pars, int model, vector<int> ages)
{
  matrix<Type> S(dim,dim);
  S.setZero();
  for(int i = 0; i< dim; i++) {
    S(i,i) = exp(2 * pars(0)); 
    if(model == 1){ //iid. nothing more to
    }
    if(model == 2){ //AR1
      for(int j = 0; j< dim; j++) {
        S(i,j) = exp(2 * pars(0)) * pow(1/(1+exp(-pars(1))), abs(ages(i)-ages(j))); //only positive correlation?
      }
    }
    if(model == 3){ //AR2: not ready
    }
  }
  return(S);
}

//this will do osa, additive/multiplicative, miss0/pool0, alternative Sigma structures (iid, AR1,...) 
template <class Type>
Type dmvnorm(vector<Type> x, vector<Type> mu, vector<Type> sig_pars, vector<int> ages, data_indicator<vector<Type>, Type> keep, int Sigma_model, 
  int do_log)
{
  matrix<Type> Sigma = make_AR_Sigma(mu.size(), sig_pars, Sigma_model, ages); //iid
  
  using namespace density;
  MVNORM_t<Type> mvnorm(Sigma);
  
  Type nll = mvnorm(x-mu, keep);
  if(do_log) return -nll;
  else return exp(-nll);

}

//usual logistic normal
//do_mult = 1: do multiplicative transformation rather than additive
template<class Type>
Type dlogisticnormal(vector<Type> obs, vector<Type> mu,  matrix<Type> S, int do_mult, int do_log)
{
  using namespace density;
  MVNORM_t<Type> mvnorm(S);

  vector<Type> x(obs.size()-1);
  if(do_mult == 1){
    x = log(obs.head(obs.size()-1));
    for(int i = 0; i < x.size(); i++) x(i) -= log(1-x.head(i+1).sum());
  }
  else x = log(obs.head(obs.size()-1)) - log(obs(obs.size()-1));
  
  Type nll = mvnorm(x-mu);
  nll += log(obs).sum(); //jacobian

  if(do_log == 1) return -nll;
  else return exp(-nll);
}

//the logistic normal with added args for osa residuals
//do_mult = 1: do multiplicative transformation rather than additive
template<class Type>
Type dlogisticnormal(vector<Type> x, vector<Type> p,  vector<Type> sig_pars, vector<int> ages, data_indicator<vector<Type>, Type> keep, int Sigma_model, 
  int do_mult, int do_log, int pool0, vector<Type> paa_obs)
{
  //NB: this is a MVN likelihood on the transformed proportions at age. 
  //MLE is equivalent, but MVN of transformed observations needed for OSA residuals.
  //x and keep should have length ages.size() - 1
  vector<Type> p_pos(ages.size());
  p_pos.setZero();
  if(pool0){ //pooling zeros
    int index = 0;
    for(int a = 0; a < p.size(); a++)
    {
      p_pos(index) += p(a);
      if((a == ages(index)-1) & (index < (ages.size() - 1))) index++;
    }
  } else { //missing zeros
    for(int a = 0; a < ages.size(); a++)
    {
      p_pos(a) = p(ages(a)-1);
    }
  }
  p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
  
  vector<Type> mu = log(p_pos.head(p_pos.size()-1));
  if(do_mult){ //multiplicative
    for(int i = 0; i < mu.size(); i++) {
      mu(i) -= log(1-p_pos.head(i+1).sum());
    }
  } else { //additive
    mu = mu - log(p_pos(p_pos.size()-1));
  }
  Type sumlog = 0;
  for(int a = 0; a < ages.size(); a++) sumlog += log(paa_obs(ages(a)-1)); //make the jacobian for the logistic normal (is this ok for OSA calculation?)
  Type ll = dmvnorm(x, mu, sig_pars, ages, keep, Sigma_model, do_log) - sumlog;
  return ll;
}


template<class Type>
Type dzinf_logisticnormal_1(vector<Type> obs, vector<Type> p, vector<Type> pars, vector<int> ages, data_indicator<vector<Type>, Type> keep, int do_log, int do_osa)
{
  //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  //probability of zero is a decreasing (logistic) function of (logit) predicted proportion caught at age.
  // 3 parameters.
  //NO OSA available!
  Type ll = 0;
  vector<Type> x = obs;
  vector<Type> X = log(p) - log(1 - p);
  //prob of zero declines with proportion caught
  vector<Type> p0 = 1.0/(1.0 + exp(exp(pars(1))*(X - pars(0)))); 
  
  int npos = 0;
  for(int a = 0; a < x.size(); a++) 
  {
    if(x(a)> 1.0e-15) {
      npos++;
      ll += log(1.0 - squeeze(p0(a))); //positives
    } else {
      ll += log(squeeze(p0(a))); //zeros
    }
  }
  if(npos>1){ //need at least two positive categories
    vector<int> pos_ind(npos);
    int k = 0;
    for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15){
      pos_ind(k) = a;
      k++;
    }
    vector<Type> x_pos = x(pos_ind);
    vector<Type> p_pos = p(pos_ind);
    p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
    //logistic normal for positive obs
    vector<Type> mu = log(p_pos.head(p_pos.size()-1)) - log(p_pos(p_pos.size()-1));
    vector<Type> Sigma_pars(1); Sigma_pars(0) = pars(2);
    matrix<Type> Sigma = make_AR_Sigma(mu.size(), Sigma_pars, 1, ages); //iid
    ll += dlogisticnormal(x_pos, mu, Sigma, 0, 1); //no osa for this one
  }
  return(ll);
}

template<class Type>
Type dzinf_logisticnormal_2(vector<Type> obs, vector<Type> p, vector<Type> pars, vector<int> ages, data_indicator<vector<Type>, Type> keep, int do_log, int do_osa)
{
  //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  //probability of zero is a decreasing function of binomial sample size with p = predicted proportion caught at age.
  // 2 parameters.
  //NO OSA available!
  Type ll = 0;
  vector<Type> x = obs;

  Type n_e = exp(pars(0)); //binomial sample size
  vector<Type> p0 = exp(n_e * log(1.0-p)); //probability of zero for the binomial distribution.
 
  int npos = 0;
  for(int a = 0; a < x.size(); a++) 
  {
    if(x(a)> 1.0e-15) {
      npos++;
      ll += log(1.0 - squeeze(p0(a))); //positives
    } else {
      ll += log(squeeze(p0(a))); //zeros
    }
  }
  if(npos>1){ //need at least two positive categories
    vector<int> pos_ind(npos);
    int k = 0;
    for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) {
      pos_ind(k) = a;
      k++;
    }
    vector<Type> x_pos = x(pos_ind);
    vector<Type> p_pos = p(pos_ind);
    p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
    //logistic normal for positive obs
    vector<Type> mu = log(p_pos.head(p_pos.size()-1)) - log(p_pos(p_pos.size()-1));
    vector<Type> Sigma_pars(1); Sigma_pars(0) = pars(1);
    matrix<Type> Sigma = make_AR_Sigma(mu.size(), Sigma_pars, 1, ages); //iid
    ll += dlogisticnormal(x_pos, mu, Sigma, 0, 1); //no osa for this one
  }
  return(ll);
}

// multivariate-Tweedie
template<class Type>
Type dmvtweedie(vector<Type> x, vector<Type> prob, Type phi, Type power, vector<int> ages, data_indicator<vector<Type>, Type> keep, int give_log=0 ){

  int n_c = x.size();
  vector<Type> p_exp(n_c);
  for(int i = 0; i < n_c; i++) p_exp(i) = prob(ages(i)-1);
  p_exp /= sum(p_exp);

  //vector<Type> p_exp(n_c);
  Type Ntotal = asDouble(sum(x)); //x.sum();
  //p_exp = prob / prob.sum();
  
  //NO OSA available
  // vector<Type> k=keep;
  // vector<int> o=order(k);
  // x=x(o); p_exp=p_exp(o); k=k(o);

  Type logres = 0;
  for( int c=0; c<n_c; c++){
    // dtweedie( Type y, Type mu, Type phi, Type p, int give_log=0 )
    logres += dtweedie( x(c), p_exp(c)*Ntotal, phi, power, true );
  }
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type get_acomp_ll(vector<Type> tf_paa_obs, vector<Type> paa_pred, Type Neff, vector<int> ages, int age_comp_model, vector<Type> age_comp_pars, 
  data_indicator<vector<Type>, Type> keep, int do_osa, vector<Type> paa_obs)
{
  Type ll = 0.0;
  vector<Type> p = paa_pred;
  if(age_comp_model == 1) {
    //multinomial
    //tf_paa_obs = Neff * paa_obs
    ll = dmultinom(tf_paa_obs, p, ages, keep, 1, 1);
  }
  if(age_comp_model == 2) {
    //saturating dirichlet-multinomial
    //tf_paa_obs = Neff * paa_obs
    vector<Type> alphas = p * exp(age_comp_pars(0));
    ll = ddirmultinom(tf_paa_obs, alphas, ages, keep, 1, 1);
  }
  if(age_comp_model == 3) { 
    //Dirichlet, miss0
    //0,1: pool 0s, do log 
    //keep, pool0, give_log, do_osa
    ll = ddirichlet(tf_paa_obs, p, exp(age_comp_pars(0)), ages, keep, 0, 1, 1);
  }
  if(age_comp_model == 4) { 
    //Dirichlet, pool0
    //0,1: pool 0s, do log 
    //keep, pool0, give_log, do_osa
    ll = ddirichlet(tf_paa_obs, p, exp(age_comp_pars(0)), ages, keep, 1, 1, 1);
  }
  if(age_comp_model == 5) { 
    //logistic-normal, miss0
    //1,0,1,0: Sigma diagonal, additive transformation, return log, missing 0s
    age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
    //need to take off obs for last age class which is NA, but keeps info on which is the last positive age
    vector<Type> x = tf_paa_obs.head(tf_paa_obs.size()-1);
    data_indicator<vector<Type>, Type> k = keep.segment(0,keep.size()-1);
    ll = dlogisticnormal(x, p, age_comp_pars, ages, k, 1, 0, 1, 0, paa_obs);
  }
  if(age_comp_model == 6) { 
    //logistic-normal, miss0, AR1 correlation
    //2,0,1,0: Sigma with AR1 cor, additive transformation, return log, missing 0s 
    age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
    vector<Type> x = tf_paa_obs.head(tf_paa_obs.size()-1);
    data_indicator<vector<Type>, Type> k = keep.segment(0,keep.size()-1);
    ll = dlogisticnormal(x, p, age_comp_pars, ages, k, 2, 0, 1, 0, paa_obs);
  }
  if(age_comp_model == 7) {
    //logistic normal. Pool zero observations with adjacent age classes.
    //1,0,1,1: Sigma diagonal, additive transformation, return log, pool 0s 
    age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
    vector<Type> x = tf_paa_obs.head(tf_paa_obs.size()-1);
    data_indicator<vector<Type>, Type> k = keep.segment(0,keep.size()-1);
    ll = dlogisticnormal(x, p, age_comp_pars, ages, k, 1, 0, 1, 1, paa_obs);
  }
  if(age_comp_model == 8) {
    //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012). 3 parameters
    //NO OSA available!
    ll = dzinf_logisticnormal_1(tf_paa_obs, p, age_comp_pars, ages, keep, 1, 0);
  }
  if(age_comp_model == 9) {
    //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
    //NO OSA available!
    ll = dzinf_logisticnormal_2(tf_paa_obs, p, age_comp_pars, ages, keep, 1, 0);
  }
  if(age_comp_model == 10) {
    //multivariate Tweedie. 2 parameters
    //NO OSA available?
    //vector<Type> temp_n = Neff * paa_obs;
    ll = dmvtweedie(tf_paa_obs, p, exp(age_comp_pars(0)), Type(1.0)+invlogit(age_comp_pars(1)), ages, keep, 1);
  }
  if(age_comp_model == 11) {
    //dirichlet-multinomial -- LINEAR
    // see https://doi.org/10.1016/j.fishres.2016.06.005
    vector<Type> alphas = sum(tf_paa_obs) * p * exp(age_comp_pars(0));
    //vector<Type> alphas = p * exp(age_comp_pars(0));
    ll = ddirmultinom(tf_paa_obs, alphas, ages, keep, 1, 1);
  }

  return ll;
}

template<class Type>
matrix<Type>  get_Neff_out(matrix<Type> Neff_input, vector<int> age_comp_models, matrix<Type> paa_pars){
  int n = age_comp_models.size();
  int n_years = Neff_input.rows();
  vector<int> any_DM(n);
  any_DM.setZero();
  for(int f = 0; f < n; f++){
    if(age_comp_models(f) == 2) any_DM(f) = 1;
    if(age_comp_models(f) == 11) any_DM(f) = 2;
  }
  matrix<Type> Neff_est = Neff_input;
  if(sum(any_DM)>0){
    //Neff_est_fleets.setZero();
    for(int f = 0; f < n; f++){
      if(any_DM(f) == 1) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) - log(N_eff) for normal D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years; y++) Neff_est(y,f) = 1 + (Neff_input(y,f) -1)/(1 + Neff_input(y,f) * exp(-paa_pars(f,0)));
      }
      if(any_DM(f) == 2) {
        // 1< N_eff_est < N_eff, logit_Neff_est = catch_paa_pars(0) for linear D-M option, so CI's could be created from that SE estimate.
        for(int y = 0; y < n_years; y++) Neff_est(y,f) = 1 + (Neff_input(y,f) -1)/(1 + exp(-paa_pars(f,0)));
      }
    }
  }
  return(Neff_est);
}
