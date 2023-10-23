template<class Type>
vector<Type> rmultinom(Type N, vector<Type> p, vector<int> ages)
{
  //multinomial
  int dim = ages.size();
  vector<Type> x(dim);
  vector<Type> p_use(dim);
  for(int i = 0; i < dim; i++) p_use(i) = p(ages(i)-1);
  p_use /= sum(p_use);
  int Nint = CppAD::Integer(N);
  x.setZero();
  for(int i = 0; i < Nint; i++)
  {
    Type y = runif(0.0,1.0);
    for(int a = 0; a < dim; a++) if(y < p_use.head(a+1).sum())
    {
      x(a) += 1.0;
      break;
    }
  }
  return x;
}

template <class Type>
vector<Type> rdirichlet(vector<Type> alpha){
  vector<Type> x=rgamma(alpha,Type(1));
  return x/sum(x);
}

template<class Type>
vector<Type> rdirmultinom(Type N, vector<Type> alpha, vector<int> ages) //dirichlet generated from iid gammas
{
  vector<Type> dp = rdirichlet(alpha);
  vector<Type> obs = rmultinom(N,dp, ages);
  return(obs);
}

template<class Type>
vector<Type> rdirichlet(vector<Type> x, vector<Type> p, Type phi, int pool0)
{
  int npos = 0;
  vector<Type> obs(x.size());
  obs.setZero();
  for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) npos++;
  vector<int> pos_ind(npos);
  vector<Type> p_pos(npos), obs_pos(npos), alpha_pos(npos);
  if(npos>1){ //need at least 2 positive categories
    int index = 0;
    for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) {
      pos_ind(index) = a;
      index++;
    }
    if(pool0){ //pooling zeros
      index = 0;
      for(int a = 0; a < x.size(); a++)
      {
        p_pos(index) += p(a);
        //x_pos(index) += x(a);
        if((x(a) > Type(1.0e-15)) & (index < npos-1)) index++;
      }
    } else { //missing zeros
      p_pos = p(pos_ind);
      //x_pos = x(pos_ind);
    }
    p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
    alpha_pos = p_pos * phi;
    obs_pos = rdirichlet(alpha_pos);
    for(int i = 0; i < npos; i++) obs(pos_ind(i)) = obs_pos(i);
  }
  else {
    obs = x;
  }
  return(obs);
}

template<class Type>
vector<Type> rdirichlet(vector<Type> p, Type phi, vector<int> ages, int pool0)
{
  vector<Type> p_pos(ages.size()), obs(ages.size());
  p_pos.setZero();
  obs.setZero();  

  if(ages.size() >1){ //need at least 2 positive categories
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
    p_pos /= sum(p_pos);//need to ensure sum(p) = 1. 

    vector<Type> alpha_pos = p_pos * phi;
    obs = rdirichlet(alpha_pos);
  }
  if(ages.size() == 1) obs(ages(0)-1) = 1.0;
  return(obs);
}

template <class Type>
vector<Type> rmvnorm(vector<Type> mu, vector<Type> sig_pars, vector<int> ages, int Sigma_model)
{
  matrix<Type> Sigma = make_AR_Sigma(mu.size(), sig_pars, Sigma_model, ages);
  
  using namespace density;
  MVNORM_t<Type> mvnorm(Sigma);
  
  vector<Type> x = mu + mvnorm.simulate();
  return(x);

}

template <class Type>
vector<Type> rmvnorm(vector<Type> p,  vector<Type> sig_pars, vector<int> ages, int Sigma_model, 
  int do_mult, int pool0)
{
  //NB: this simulates a MVN rv with mean defined by logit-normal proportions at age. 
  //the simulated vector will be MVN with length p.size() - 1
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
  p_pos /= sum(p_pos); //necessary for multiplicitive
    
  vector<Type> mu = log(p_pos.head(p_pos.size()-1));
  if(do_mult){ //multiplicative
    for(int i = 0; i < mu.size(); i++) {
      mu(i) -= log(1-p_pos.head(i+1).sum());
    }
  } else { //additive
    mu = mu - log(p_pos(p_pos.size()-1));
  }

  vector<Type> x = rmvnorm(mu, sig_pars, ages, Sigma_model);
  return x;
}

template <class Type>
vector<Type> mvn_to_LN(vector<Type> x, int do_mult)
{
  int d = x.size();
  vector<Type> y(d+1);

  if(do_mult){
    Type denom = 1;
    for(int i = 0; i < d; i++) {
      denom *= 1 + exp(x(i));
      y(i) = exp(x(i))/denom;
    }
    y(d) = 1/denom;
  } else {
    for(int i = 0; i < d; i++) y(i) = exp(x(i))/(1 + sum(exp(x)));
    y(d) = 1/(1 + sum(exp(x)));
  }
  return(y);
}

template <class Type>
vector<Type> rlogisticnormal(vector<Type> x, vector<Type> p, vector<Type> pars, vector<int> ages, int Sigma_model, int do_mult, int pool0)
{
  int npos = 0;
  vector<Type> obs(x.size());
  obs.setZero();
  for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) npos++;
  vector<int> pos_ind(npos);
  vector<Type> p_pos(npos);
  matrix<Type> Sigma_full = make_AR_Sigma(p.size(), pars, Sigma_model, ages); //iid

  if(npos>1){ //need at least 2 positive categories
    int index = 0;
    for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) {
      pos_ind(index) = a;
      index++;
    }
    matrix<Type> Sigma(npos-1, npos-1); //reduced dimensions for npos-1 values
    for(int i = 0; i < npos - 1; i++){
      for(int j = 0; j < npos - 1; j++){
        Sigma(i,j) = Sigma_full(pos_ind(i),pos_ind(j));
      }
    }

    if(pool0){ //pooling zeros
      index = 0;
      for(int a = 0; a < x.size(); a++)
      {
        p_pos(index) += p(a);
        if((x(a) > Type(1.0e-15)) & (index < npos-1)) index++;
      }
    } else { //missing zeros
      p_pos = p(pos_ind);
    }
    p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
    vector<Type> mu = log(p_pos.head(p_pos.size()-1));
    if(do_mult){ //multiplicative
      for(int i = 0; i < mu.size(); i++) {
        mu(i) -= log(1-p_pos.head(i+1).sum());
      }
    }else { //additive
      mu = mu - log(p_pos(p_pos.size()-1));
    }

    using namespace density;
    MVNORM_t<Type> mvnorm(Sigma);
    vector<Type> y = exp(mu + mvnorm.simulate());
    if(do_mult){
      for(int i = 0; i < npos-1; i++) y = y/(1 + (y.head(i+1))).prod();
    } else {
      y = y/(1 + sum(y));
    }
    for(int i = 0; i < npos-1; i++) obs(pos_ind(i)) = y(i);
    obs(pos_ind(npos-1)) = 1 - sum(y);
  }
  else {
    obs = x;
  }
  return(obs);  
}
template<class Type>
vector<Type> rzinf_logisticnormal_1(vector<Type> p, vector<Type> pars, vector<int> ages)
//Type zinf_logisticnormal_1(vector<Type> obs, vector<Type> p, vector<Type> pars, vector<Type> keep, vector<Type> keep_l, vector<Type> keep_h, int do_log, int do_osa)
{
  //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  //probability of zero is a decreasing (logistic) function of (logit) predicted proportion caught at age.
  // 3 parameters.
  //NO OSA available!
  //Type ll = 0;
  vector<Type> obs(p.size());
  obs.setZero();
  vector<Type> X = log(p) - log(1 - p);
  //prob of zero declines with proportion caught
  vector<Type> p0 = 1.0/(1.0 + exp(exp(pars(1))*(X - pars(0)))); 
  for(int a = 0; a < obs.size(); a++) obs(a) = rbinom(Type(1.0), Type(1.0) - p0(a)); // generate instances of positive observations
  
  //vector<Type> x = obs;
  int npos = 0;
  for(int a = 0; a < obs.size(); a++) if(obs(a) > 1.0e-15) npos++;
  if(npos>1) //need at least 2 positive obs
  {
    vector<Type> Sigma_pars(1); Sigma_pars(0) = pars(2);
    obs = rlogisticnormal(obs, p, Sigma_pars, ages, 1, 0, 0); //1,0,0: 
    //obs = obs_out;
  }
  return(obs);
}

template<class Type>
vector<Type> rzinf_logisticnormal_2(vector<Type> p, vector<Type> pars, vector<int> ages)
{
  //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  //probability of zero is a decreasing function of binomial sample size with p = predicted proportion caught at age.
  // 2 parameters.
  vector<Type> obs(p.size());
  obs.setZero();
  Type n_e = exp(pars(0)); //binomial sample size
  vector<Type> p0 = exp(n_e * log(1.0-p)); //probability of zero for the binomial distribution.

  for(int a = 0; a < obs.size(); a++) obs(a) = rbinom(Type(1.0), Type(1.0) - p0(a)); // generate instances of positive observations
  
  //vector<Type> x = obs;
  int npos = 0;
  for(int a = 0; a < obs.size(); a++) if(obs(a) > 1.0e-15) npos++;
  if(npos>1) //need at least 2 positive obs
  {
    vector<Type> Sigma_pars(1); Sigma_pars(0) = pars(1);
    obs = rlogisticnormal(obs, p, Sigma_pars, ages, 1, 0, 0); //Sigma_model =iid, additive transform, miss0
    //obs = obs_out;
  }
  return(obs);
}

// Simulate from tweedie
// Adapted from tweedie::rtweedie function in R
template<class Type>
Type rTweedie( Type mu, Type phi, Type power){
 Type lambda = pow(mu, Type(2.0) - power) / (phi * (Type(2.0) - power));
 Type alpha = (Type(2.0) - power) / (Type(1.0) - power);
 Type gam = phi * (power - Type(1.0)) * pow(mu, power - Type(1.0));
 Type N = rpois(lambda);
 Type B = rgamma(-N * alpha, gam);   /// Using Shape-Scale parameterization
 return B;
}

template<class Type>
vector<Type> rmvtweedie( Type N, vector<Type> p, Type phi, Type power, vector<int> ages)
{
  int dim = ages.size();
  vector<Type> obs(dim);
  //obs.setZero();
  for(int i = 0; i < dim; i++)
  {
    obs(i) = rTweedie( N*p(ages(i)-1), phi, power );
  }
  return(obs);
}

template<class Type>
vector<Type> sim_acomp(vector<Type> paa_pred, Type Neff, vector<int> ages, int age_comp_model, vector<Type> age_comp_pars)
{
  int n_ages = ages.size();
  vector<Type> obs(n_ages);
  vector<Type> p = paa_pred;
  obs.setZero();
  if(age_comp_model == 1)
  {
    //multinomial
    p += 1.0e-15; //for log of any p = 0
    obs = rmultinom(Neff, p, ages);
    //obs = obs/obs.sum();// proportions
  }
  if(age_comp_model == 2) //dirichlet-multinomial. dirichlet generated from iid gammas and multinomial from uniform
  {
    //int N = CppAD::Integer(Neff);
    vector<Type> alpha = p * exp(age_comp_pars(0));
    obs = rdirmultinom(Neff,alpha,ages);
    //obs = obs/obs.sum();// proportions
  }
  if(age_comp_model == 3) { //Dirichlet, miss0
    obs = rdirichlet(p, exp(age_comp_pars(0)), ages, 0); //0: miss0
  }
  if(age_comp_model == 4) { //Dirichlet, pool0
    obs = rdirichlet(p, exp(age_comp_pars(0)), ages, 1); //1: pool0
  }
  if(age_comp_model == 5) { //logistic-normal, miss0
    age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
    obs = rmvnorm(p, age_comp_pars, ages, 1, 0, 0); //1,0,0: Sigma diagonal, additive transformation, missing 0s 
  }
  if(age_comp_model == 6) { //logistic-normal, ar1 cor, miss0
    age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
    obs = rmvnorm(p, age_comp_pars, ages, 2, 0, 0); //2,0,0: Sigma AR1 cor, additive transformation, missing 0s 
  }
  if(age_comp_model == 7) { //logistic-normal, pool0
    age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
    obs = rmvnorm(p, age_comp_pars, ages, 1, 0, 1); //1,0,1: Sigma diagonal, additive transformation, pooling 0s 
  }
  if(age_comp_model == 8) 
  {
    //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
    //NO OSA available!
    obs = rzinf_logisticnormal_1(p, age_comp_pars, ages);
  }
  if(age_comp_model == 9) 
  {
    //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
    //NO OSA available!
    obs = rzinf_logisticnormal_2(p, age_comp_pars, ages);
  }
  if(age_comp_model == 10) 
  {
    //multivariate Tweedie. 2 parameters
    obs = rmvtweedie(Neff, p, exp(age_comp_pars(0)), Type(1.0)+invlogit(age_comp_pars(1)), ages);
  }
  if(age_comp_model == 11) //"linearized" dirichlet-multinomial. dirichlet generated from iid gammas and multinomial from uniform
  {
    //int N = CppAD::Integer(Neff);
    vector<Type> alpha = Neff * p * exp(age_comp_pars(0));
    obs = rdirmultinom(Neff,alpha,ages);
    //obs = obs/obs.sum();// proportions
  }
  return obs;
}

//make proporportions at age observations from transformed versions
template<class Type>
vector<Type> make_paa(vector<Type> tf_paa_obs, int age_comp_model, vector<int> ages, int n_ages)
{
  //int n_ages = paa_obs.size();
  vector<Type> paa_out(n_ages);
  paa_out.setZero();
  if((age_comp_model <5) | (age_comp_model > 7)) {
    for(int i = 0; i < ages.size(); i++) paa_out(ages(i)-1) = tf_paa_obs(i); //identity transform, zeros allowed
    if((age_comp_model < 3) | (age_comp_model > 9)) paa_out /= sum(paa_out); //multinomial, D-M, mvtweedie are in numbers, linear D-M
  }
  if((age_comp_model > 4) & (age_comp_model < 8)) { //logistic normal
    vector<Type> p_pos(ages.size());
    if(tf_paa_obs.size() == ages.size()){ //unused last element has been provided
      p_pos = mvn_to_LN(vector<Type> (tf_paa_obs.head(tf_paa_obs.size()-1)), 0);// no multiplicative options right now
    } else {
      p_pos = mvn_to_LN(tf_paa_obs, 0);// no multiplicative options right now
    }
    for(int i = 0; i < ages.size(); i++) paa_out(ages(i)-1) = p_pos(i); 
  }
  return paa_out;
}
