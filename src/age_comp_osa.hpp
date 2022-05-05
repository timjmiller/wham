//from Nielsen et al. (in prep.). OSA residuals for composition data
//osa residual versions of dmultinom, ddirmultinom, dlogisticnormal
namespace my_atomic {
  /*
   *  Modified from R source pbinom.c
   *  Mathlib : A C Library of Special Functions
   *  Copyright (C) 1998 Ross Ihaka
   *  Copyright (C) 2000-2015  The R Core Team
   *  Copyright (C) 2004-2015  The R Foundation
   */
  template<class T> int R_finite(T x) { return std::isfinite(asDouble(x)); }
  template<class T> int isnan(T x) { return std::isnan(asDouble(x)); }

#undef ML_ERROR
#undef MATHLIB_ERROR
#undef MATHLIB_WARNING
#undef MATHLIB_WARNING2
#undef MATHLIB_WARNING3
#undef MATHLIB_WARNING4
#undef MATHLIB_WARNING5
#undef ML_POSINF
#undef ML_NEGINF
#undef ML_NAN
#undef M_SQRT_2dPI
#undef ISNAN
# define ML_ERROR(x, s) /* nothing */
# define MATHLIB_ERROR(fmt,x) /* nothing */
# define MATHLIB_WARNING(fmt,x) /* nothing */
# define MATHLIB_WARNING2(fmt,x,x2) /* nothing */
# define MATHLIB_WARNING3(fmt,x,x2,x3) /* nothing */
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) /* nothing */
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) /* nothing */
#define ML_POSINF R_PosInf
#define ML_NEGINF R_NegInf
#define ML_NAN    R_NaN
#define M_SQRT_2dPI 0.797884560802865355879892119869  /* sqrt(2/pi) */
#define ISNAN(x) (isnan(x)!=0)

#define ML_ERR_return_NAN return R_NaN

#define R_D__0  (log_p ? ML_NEGINF : 0.)    /* 0 */
#define R_D__1  (log_p ? 0. : 1.)     /* 1 */
#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)    /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0) /* 1 */
# define attribute_hidden __attribute__ ((visibility ("hidden")))


  template<class Float>
  attribute_hidden
  Float pbinom0_raw(Float x, Float n, Float p, Float lower_tail_, Float log_p_)
  {
    int lower_tail = (int)trunc(lower_tail_);
    int log_p = (int)trunc(log_p_);

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n) || ISNAN(p))
      return x + n + p;
    if (!R_FINITE(n) || !R_FINITE(p)) ML_ERR_return_NAN;

#endif
    // if(R_nonint(n)) {
    //   MATHLIB_WARNING(("non-integer n = %f"), n);
    //   //ML_ERR_return_NAN;
    //   return 0;
    // }
    n = (Float)(int)trunc(n);
    /* PR#8560: n=0 is a valid value */
    if(n < 0 || p < 0 || p > 1) ML_ERR_return_NAN;

    if (x < 0) return R_DT_0;
    //x = floor(x + 1e-7);
    if (n <= x) return R_DT_1;
    return atomic::toms708::pbeta((Float)p, (Float)(x + 1), (Float)(n - x), (int)!lower_tail, (int)log_p);
  }
  
  template<class Float>
  Float pbinom0(Float x, Float n, Float p, Float lower_tail, Float log_p)
  {
    return pbinom0_raw(x,n,p,lower_tail,log_p);
  }

  TMB_BIND_ATOMIC(pbinom1,11100,pbinom0(x[0], x[1], x[2], x[3], x[4]))

}

template<class Type>
Type pbinom(Type x, Type n, Type p, int lower_tail, int log_p){
  CppAD::vector<Type> tx(6);
  tx[0] = x;
  tx[1] = n;
  tx[2] = p;
  tx[3] = lower_tail;
  tx[4] = log_p;
  tx[5] = 0; // extra argument for derivative order
  Type res = my_atomic::pbinom1(tx)[0];
  return res;
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
Type dmultinom(vector<Type> x, vector<Type> p, vector<Type> keep, vector<Type> keep_l, vector<Type> keep_h, int give_log, int do_osa)
{
  Type logres = 0;
  if(do_osa){
    vector<Type> k=keep;
    vector<Type> l=keep_l;
    vector<Type> h=keep_h;
    vector<int> o=order(k);
    vector<Type> kno = k;
    x=x(o); p=p(o); k=k(o); l=l(o); h=h(o);
    Type nUnused=sum(x);
    Type pUsed=0;
    Type cdf;
    for(int i=0; i<x.size(); ++i){
      if(i!=(x.size()-1)){
        logres += k(i)*dbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true);
        cdf = pbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true,false);
        nUnused -= x(i);
        pUsed += p(i);
      }else{ // last index 
        logres += k(i)*Type(0);
        cdf = Type(1);
      }
      cdf = squeeze(cdf);
      logres += l[i] * log( cdf );       // NaN protected
      logres += h[i] * log( 1.0 - cdf ); // NaN protected
    }
  } else {
    logres = dmultinom(x,p,1);
  }
  if(give_log){
    return logres;
  }else{ 
    return exp(logres);
  }
}

template <class Type>
Type dmultinom(vector<Type> x, vector<Type> p, data_indicator<vector<Type>, Type> keep, int give_log, int do_osa)
{
  Type logres = 0;
  if(do_osa){
    vector<Type> k=keep;
    vector<Type> l=keep.cdf_lower;
    vector<Type> h=keep.cdf_upper;
    vector<int> o=order(k);
    vector<Type> kno = k;
    x=x(o); p=p(o); k=k(o); l=l(o); h=h(o);
    Type nUnused=sum(x);
    Type pUsed=0;
    Type cdf;
    for(int i=0; i<x.size(); ++i){
      if(i!=(x.size()-1)){
        logres += k(i)*dbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true);
        cdf = pbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true,false);
        nUnused -= x(i);
        pUsed += p(i);
      }else{ // last index 
        logres += k(i)*Type(0);
        cdf = Type(1);
      }
      cdf = squeeze(cdf);
      logres += l[i] * log( cdf );       // NaN protected
      logres += h[i] * log( 1.0 - cdf ); // NaN protected
    }
  } else {
    logres = dmultinom(x,p,1);
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
Type ddirichlet(vector<Type> x, vector<Type> alpha, data_indicator<vector<Type>, Type> keep, int give_log, int do_osa)
{
  Type logres = 0;
  if(do_osa){
    vector<Type> k=keep;
    vector<Type> l=keep.cdf_lower;
    vector<Type> h=keep.cdf_upper;
    vector<int> o=order(k);
    x=x(o); alpha=alpha(o); k=k(o); l=l(o); h=h(o);  
    int n = alpha.size();
    Type cdf;
    Type sx = x.sum();
    Type sa = alpha.sum();
    sa -= alpha(0);
    logres = k(0)*dbeta(x(0),alpha(0),sa,true);
    cdf = pbeta(x(0),alpha(0),sa);
    cdf = squeeze(cdf);
    logres += l(0) * log( cdf );       
    logres += h(0) * log( 1.0 - cdf ); 
    
    for(int i=1; i<(n-1); ++i){
      sx -= x(i-1);
      sa -= alpha(i);
      logres += k(i)*(dbeta(x(i)/sx,alpha(i),sa,true)-log(sx));
      cdf = pbeta(x(i)/sx,alpha(i),sa);
      cdf = squeeze(cdf);
      logres += l(i) * log( cdf );       
      logres += h(i) * log( 1.0 - cdf ); 
    }
    logres += k(n-1)*Type(0);
    cdf=Type(1);
    cdf = squeeze(cdf);
    logres += l(n-1) * log( cdf );       
    logres += h(n-1) * log( 1.0 - cdf ); 
  }
  else{
    logres = ddirichlet(x,alpha,1);
  }
  
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}

template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> p, Type phi, data_indicator<vector<Type>, Type> keep, int pool0, int give_log, int do_osa)
{
  Type logres = 0;
  int npos = 0;
  for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) npos++;
  vector<int> pos_ind(npos);
  vector<Type> p_pos(npos), x_pos(npos), alpha_pos(npos);
  Type cdf;

  if(npos>1){ //need at least 2 positive categories
    if(do_osa){
      vector<Type> k=keep;
      vector<Type> l=keep.cdf_lower;
      vector<Type> h=keep.cdf_upper;
      vector<int> o=order(k);
      x=x(o); p=p(o); k=k(o); l=l(o); h=h(o);  

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
          x_pos(index) += x(a);
          if(x(a) > Type(1.0e-15) & index < npos-1) index++;
        }
      } else { //missing zeros
        p_pos = p(pos_ind);
        x_pos = x(pos_ind);
      }
      alpha_pos = p_pos * phi;

      int n = alpha_pos.size();
      Type sx = x_pos.sum();
      Type sa = alpha_pos.sum();
      sa -= alpha_pos(0);
      logres += k(pos_ind(0))*dbeta(x_pos(0),alpha_pos(0),sa,true);
      cdf = squeeze(pbeta(x_pos(0),alpha_pos(0),sa));
      logres += l(pos_ind(0)) * log( cdf );       
      logres += h(pos_ind(0)) * log( 1.0 - cdf ); 
      for(int i=1; i<(n-1); ++i){
        sx -= x_pos(i-1);
        sa -= alpha_pos(i);
        logres += k(pos_ind(i))*(dbeta(x_pos(i)/sx,alpha_pos(i),sa,true)-log(sx));
        cdf = squeeze(pbeta(x_pos(i)/sx,alpha_pos(i),sa));
        logres += l(pos_ind(i)) * log( cdf );       
        logres += h(pos_ind(i)) * log( 1.0 - cdf ); 
      }
      //fill in the zeros and last positive
      cdf = squeeze(Type(1));
      logres += k(pos_ind(npos-1)) * Type(0);
      logres += l(pos_ind(npos-1)) * log( cdf );       
      logres += h(pos_ind(npos-1)) * log( 1.0 - cdf ); 
      for(int i=1; i< x.size(); ++i) if(x(i) < 1.0e-15) {
        logres += k(i) * Type(0);
        logres += l(i) * log( cdf );       
        logres += h(i) * log( 1.0 - cdf ); 
      }

    }
    else{ //do likelihood without osa machinery
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
          x_pos(index) += x(a);
          if(x(a) > Type(1.0e-15) & index < npos-1) index++;
        }
      } else { //missing zeros
        p_pos = p(pos_ind);
        x_pos = x(pos_ind);
      }
      p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
      alpha_pos = p_pos * phi;
      logres = ddirichlet(x_pos,alpha_pos,1);
    }
  } 
  else {
    cdf = squeeze(Type(1));
    for(int i = 0; i < keep.size(); i++) {
      logres += keep(i) * Type(0); 
      logres += keep.cdf_lower(i) * log(cdf);
      logres += keep.cdf_upper(i) * log(1.0 - cdf); 
    }
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
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, data_indicator<vector<Type>, Type> keep, int do_log, int do_osa)
//Type ddirmultinom(vector<Type> obs, vector<Type> alpha, vector<Type> keep, vector<Type> keep_l, vector<Type> keep_h, int do_log, int do_osa)
{
  int dim = obs.size();
  Type N = obs.sum(), ll = 0.0;
  if(do_osa){
    vector<Type> k=keep;
    vector<Type> l=keep.cdf_lower;
    vector<Type> h=keep.cdf_upper;
    vector<int> o=order(k);
    obs=obs(o); alpha=alpha(o); k=k(o); l=l(o); h=h(o);   
    vector<Type> alphas_a(2), obs_a(2);
    for(int a = 1; a < dim; a++){
      obs_a(0) = obs(a-1);
      obs_a(1) = obs.tail(dim-a).sum();
      alphas_a(0) = alpha(a-1);
      alphas_a(1) = alpha.tail(dim-a).sum();
      ll += k(a-1) * ddirmultinom(obs_a, alphas_a, 1); //beta-binomial, just two categories
      Type cdf = pbetabinom(obs_a(0), obs_a.sum(), alphas_a(0), alphas_a(1),0);
      cdf = squeeze(cdf);
      ll += l(a-1) * log( cdf );       
      ll += h(a-1) * log( 1.0 - cdf ); 
    }
    ll += k(dim-1)*Type(0);
    Type cdf= squeeze(Type(1));
    ll += l(dim-1) * log( cdf );       
    ll += h(dim-1) * log( 1.0 - cdf ); 
  }
  else{
    ll = ddirmultinom(obs,alpha,1);
  }
  if(do_log == 1) return ll;
  else return exp(ll);
}

//the logistic normal with added args for osa residuals
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


//function to generate covariance matrix for logistic normal
template<class Type>
matrix<Type> make_LN_Sigma(int dim, vector<Type> pars, int model){
  matrix<Type> S(dim,dim);
  S.setZero();
  for(int i = 0; i< dim; i++) {
    S(i,i) = exp(2 * pars(0)); 
    if(model == 1){ //iid. nothing more to
    }
    if(model == 2){ //AR1: not ready
      for(int j = 0; j< dim; j++) {
        S(i,j) = exp(2 * pars(0)) * pow(1/(1+exp(-pars(1))), abs(i-j)); //only positive correlation?
      }
    }
    if(model == 3){ //AR2: not ready
    }
  }
  return(S);
}

//this will do osa, additive/multiplicative, miss0/pool0, alternative Sigma structures (iid, AR1,...) 
template <class Type>
Type dlogisticnormal(vector<Type> x, vector<Type> p, vector<Type> pars, data_indicator<vector<Type>, Type> keep, int Sigma_model, 
  int do_mult, int give_log, int pool0, int do_osa)
{
  Type logres = 0;
  int npos = 0;
  for(int a = 0; a < x.size(); a++) if(x(a)> 1.0e-15) npos++;
  vector<Type> x_pos(npos), p_pos(npos);
  x_pos.setZero(); p_pos.setZero();
  vector<int> pos_ind(npos);
  Type cdf;
  matrix<Type> Sigma_full = make_LN_Sigma(p.size(), pars, Sigma_model); //iid

  if(npos>1){ //need at least 2 positive categories
    if(do_osa){
      vector<Type> k=keep;
      vector<Type> l=keep.cdf_lower;
      vector<Type> h=keep.cdf_upper;
      vector<int> o=order(k);
      x=x(o); p=p(o); k=k(o); l=l(o); h=h(o);  

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
          x_pos(index) += x(a);
          if(x(a) > Type(1.0e-15) & index < npos-1) index++;
        }
      } else { //missing zeros
        p_pos = p(pos_ind);
        x_pos = x(pos_ind);
      }

      p_pos = p_pos/p_pos.sum(); //p rescaled to sum to 1 for observed age classes
      //additive transformation
      vector<Type> mu = log(p_pos.head(p_pos.size()-1));// - log(p_pos(p_pos.size()-1));
      vector<Type> y = log(x_pos.head(x_pos.size()-1));
      if(do_mult){ //multiplicative
        for(int i = 0; i < y.size(); i++) {
          y(i) -= log(1-x_pos.head(i+1).sum());
          mu(i) -= log(1-p_pos.head(i+1).sum());
        }
      }else { //additive
        y = y - log(x_pos(x_pos.size()-1));
        mu = mu - log(p_pos(p_pos.size()-1));
      }
      
      //first positive value
      Type m = mu(0);
      Type sd = pow(Sigma(0,0),0.5);
      logres += k(pos_ind(0))*dnorm(y(0),m, sd,true);
      cdf = squeeze(pnorm(y(0),m,sd));
      logres += l(pos_ind(0)) * log( cdf );       
      logres += h(pos_ind(0)) * log( 1.0 - cdf ); 
      
      //all the information is in the n-1 positive classes.
      for(int i=1; i<(npos-1); ++i){ //only positive values
        matrix<Type> S_other(i,i); 
        matrix<Type> S_row(1,i); 
        matrix<Type> S_col(i,1);
        matrix<Type> res_other(i,1);
        for(int j = 0; j < i; j++) {
          S_row(0,j) = Sigma(i,j);
          S_col(j,0) = Sigma(j,i);
          res_other(j,0) = y(j) - mu(j);
          for(int a = 0; a < i; a++) S_other(j,a) = Sigma(j,a);
        }
        matrix<Type> inv_S_other = atomic::matinv(S_other);
        m = mu(i) - (S_row * inv_S_other * res_other)(0,0);
        //mu_cond = mu_i - S_i,other %*% S^-1 %*% (y_other - mu_other)
        sd = Sigma(i,i) - (S_row * inv_S_other * S_col)(0,0);
        sd = pow(sd, 0.5);
        logres += k(pos_ind(i)) * dnorm(y(i),m,sd,true);
        cdf = squeeze(pnorm(y(i), m, sd));
        logres += l(pos_ind(i)) * log( cdf );       
        logres += h(pos_ind(i)) * log( 1.0 - cdf ); 
      }
      logres += k(pos_ind(npos-1)) * Type(0);
      cdf = squeeze(Type(1));
      logres += l(pos_ind(npos-1)) * log( cdf );       
      logres += h(pos_ind(npos-1)) * log( 1.0 - cdf ); 
      
      //dealing with the zero classes
      for(int i = 0; i < x.size(); i++) if(x(i) < 1.0e-15) {
        logres += k(i)*Type(0); 
        cdf = squeeze(Type(1));
        logres += l(i) * log(cdf);
        logres += h(i) * log(1.0 - cdf); 
      }
    }
    else{ //do likelihood without osa machinery
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
          x_pos(index) += x(a);
          if(x(a) > Type(1.0e-15) & index < npos-1) index++;
        }
      } else { //missing zeros
        p_pos = p(pos_ind);
        x_pos = x(pos_ind);
      }
      p_pos /= p_pos.sum(); //need to ensure sum(p) = 1. x should already sum to 1.
      vector<Type> mu = log(p_pos.head(p_pos.size()-1));// - log(p_pos(p_pos.size()-1));
      vector<Type> y = log(x_pos.head(x_pos.size()-1));
      if(do_mult){ //multiplicative
        for(int i = 0; i < y.size(); i++) {
          y(i) -= log(1-x_pos.head(i+1).sum());
          mu(i) -= log(1-p_pos.head(i+1).sum());
        }
      }else { //additive
        y = y - log(x_pos(x_pos.size()-1));
        mu = mu - log(p_pos(p_pos.size()-1));
      }

      using namespace density;
      MVNORM_t<Type> mvnorm(Sigma);
      logres += -mvnorm(y-mu);
      logres -= log(x_pos).sum(); //jacobian, do it only when osa residuals not being calculated.
    }
  }
  else {
    cdf = squeeze(Type(1));
    for(int i = 0; i < keep.size(); i++) {
      logres += keep(i) * Type(0); 
      logres += keep.cdf_lower(i) * log(cdf);
      logres += keep.cdf_upper(i) * log(1.0 - cdf); 
    }
  }
  
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}

template<class Type>
Type dzinf_logisticnormal_1(vector<Type> obs, vector<Type> p, vector<Type> pars, data_indicator<vector<Type>, Type> keep, int do_log, int do_osa)
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
    matrix<Type> Sigma = make_LN_Sigma(mu.size(), Sigma_pars, 1); //iid
    ll += dlogisticnormal(x_pos, mu, Sigma, 0, 1); //no osa for this one
  }
  return(ll);
}

template<class Type>
Type dzinf_logisticnormal_2(vector<Type> obs, vector<Type> p, vector<Type> pars, data_indicator<vector<Type>, Type> keep, int do_log, int do_osa)
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
    matrix<Type> Sigma = make_LN_Sigma(mu.size(), Sigma_pars, 1); //iid
    ll += dlogisticnormal(x_pos, mu, Sigma, 0, 1); //no osa for this one
  }
  return(ll);
}

template<class Type>
Type get_acomp_ll(vector<Type> paa_obs, vector<Type> paa_pred, Type Neff, int age_comp_model, vector<Type> age_comp_pars, 
  data_indicator<vector<Type>, Type> keep, 
  int do_osa)
{
  Type ll = 0.0;
  int n_ages = paa_obs.size();
  vector<Type> p = paa_pred + 1.0e-15;
  //if(use_obs){
    if(age_comp_model == 1) //multinomial
    {
      vector<Type> x = Neff * paa_obs;
      ll = dmultinom(x, p, keep, 1, do_osa);
    }
    if(age_comp_model == 2) //dirichlet-multinomial
    {
      vector<Type> x = Neff * paa_obs;
      vector<Type> alphas = p * exp(age_comp_pars(0));
      ll = ddirmultinom(x, alphas, keep, 1, do_osa);
    }
    if(age_comp_model == 3) { //Dirichlet, miss0
      //0,1: pool 0s, do log 
      //keep, pool0, give_log, do_osa
      ll = ddirichlet(paa_obs, p, exp(age_comp_pars(0)), keep, 0,1,do_osa);
    }
    if(age_comp_model == 4) { //Dirichlet, pool0
      //0,1: pool 0s, do log 
      //keep, pool0, give_log, do_osa
      ll = ddirichlet(paa_obs, p, exp(age_comp_pars(0)), keep, 1,1,do_osa);
    }
    if(age_comp_model == 5) { //logistic-normal, miss0
      //1,0,1,0: Sigma diagonal, additive transformation, return log, missing 0s
      age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
      ll = dlogisticnormal(paa_obs, p, age_comp_pars, keep, 1, 0, 1, 0, do_osa);
    }
    if(age_comp_model == 6) { //logistic-normal, miss0, AR1 correlation
      //2,0,1,0: Sigma with AR1 cor, additive transformation, return log, missing 0s 
      age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
      ll = dlogisticnormal(paa_obs, p, age_comp_pars, keep, 2, 0, 1, 0, do_osa);
    }
    if(age_comp_model == 7) 
    {
      //logistic normal. Pool zero observations with adjacent age classes.
      //1,0,1,1: Sigma diagonal, additive transformation, return log, pool 0s 
      age_comp_pars(0) -= 0.5*log(Neff); //an adjustment for interannual variation in sampling effort
      ll = dlogisticnormal(paa_obs, p, age_comp_pars, keep, 1, 0, 1, 1, do_osa);
    }
    if(age_comp_model == 8) 
    {
      //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012). 3 parameters
      //NO OSA available!
      ll = dzinf_logisticnormal_1(paa_obs, p, age_comp_pars, keep, 1, 0);
    }
    if(age_comp_model == 9) 
    {
      //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
      //NO OSA available!
      ll = dzinf_logisticnormal_2(paa_obs, p, age_comp_pars, keep, 1, 0);
    }

  return ll;
}

template<class Type>
Type get_lcomp_ll(vector<Type> pal_obs, vector<Type> pal_pred, Type NeffL, int len_comp_model, vector<Type> len_comp_pars, 
  data_indicator<vector<Type>, Type> keep, 
  int do_osa)
{
  Type ll = 0.0;
  int n_lengths = pal_obs.size();
  vector<Type> p = pal_pred + 1.0e-15;
  //if(use_obs){
    if(len_comp_model == 1) //multinomial
    {
      vector<Type> x = NeffL * pal_obs;
      ll = dmultinom(x, p, keep, 1, do_osa);
    }
    if(len_comp_model == 2) //dirichlet-multinomial
    {
      vector<Type> x = NeffL * pal_obs;
      vector<Type> alphas = p * exp(len_comp_pars(0));
      ll = ddirmultinom(x, alphas, keep, 1, do_osa);
    }
    if(len_comp_model == 3) { //Dirichlet, miss0
      //0,1: pool 0s, do log 
      //keep, pool0, give_log, do_osa
      ll = ddirichlet(pal_obs, p, exp(len_comp_pars(0)), keep, 0,1,do_osa);
    }
    if(len_comp_model == 4) { //Dirichlet, pool0
      //0,1: pool 0s, do log 
      //keep, pool0, give_log, do_osa
      ll = ddirichlet(pal_obs, p, exp(len_comp_pars(0)), keep, 1,1,do_osa);
    }
    if(len_comp_model == 5) { //logistic-normal, miss0
      //1,0,1,0: Sigma diagonal, additive transformation, return log, missing 0s
      len_comp_pars(0) -= 0.5*log(NeffL); //an adjustment for interannual variation in sampling effort
      ll = dlogisticnormal(pal_obs, p, len_comp_pars, keep, 1, 0, 1, 0, do_osa);
    }
    if(len_comp_model == 6) { //logistic-normal, miss0, AR1 correlation
      //2,0,1,0: Sigma with AR1 cor, additive transformation, return log, missing 0s 
      len_comp_pars(0) -= 0.5*log(NeffL); //an adjustment for interannual variation in sampling effort
      ll = dlogisticnormal(pal_obs, p, len_comp_pars, keep, 2, 0, 1, 0, do_osa);
    }
    if(len_comp_model == 7) 
    {
      //logistic normal. Pool zero observations with adjacent age classes.
      //1,0,1,1: Sigma diagonal, additive transformation, return log, pool 0s 
      len_comp_pars(0) -= 0.5*log(NeffL); //an adjustment for interannual variation in sampling effort
      ll = dlogisticnormal(pal_obs, p, len_comp_pars, keep, 1, 0, 1, 1, do_osa);
    }
    if(len_comp_model == 8) 
    {
      //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012). 3 parameters
      //NO OSA available!
      ll = dzinf_logisticnormal_1(pal_obs, p, len_comp_pars, keep, 1, 0);
    }
    if(len_comp_model == 9) 
    {
      //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
      //NO OSA available!
      ll = dzinf_logisticnormal_2(pal_obs, p, len_comp_pars, keep, 1, 0);
    }

  return ll;
}