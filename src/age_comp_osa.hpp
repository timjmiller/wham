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
Type dmultinom_osa(vector<Type> x, vector<Type> p, data_indicator<vector<Type>, Type> keep, int give_log)
{
  vector<Type> k=keep;
  vector<Type> l=keep.cdf_lower;
  vector<Type> h=keep.cdf_upper;
  vector<int> o=order(k);
  x=x(o); p=p(o); k=k(o); l=l(o); h=h(o);
  Type logres=0;
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
  if(give_log){
    return logres;
  }else{ 
    return exp(logres);
  }
}

template <class Type>
Type ddirichlet_osa(vector<Type> x, vector<Type> p, Type phi, data_indicator<vector<Type>, Type> keep, int give_log)
{
  vector<Type> x = obs;
  vector<Type> alpha = (p + Type(10e-15)) * phi;
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
  Type logres=k(0)*dbeta(x(0),alpha(0),sa,true);
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
  
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
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
Type ddirmultinom_osa(vector<Type> obs, vector<Type> alpha, data_indicator<vector<Type>, Type> keep, int do_log = 0)
{

  vector<Type> k=keep;
  vector<Type> l=keep.cdf_lower;
  vector<Type> h=keep.cdf_upper;

  vector<int> o=order(k);
  obs=obs(o); alpha=alpha(o); k=k(o); l=l(o); h=h(o);
 
  int dim = obs.size();
  Type N = obs.sum(), ll = 0.0;
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
  if(do_log == 1) return ll;
  else return exp(ll);
}

//the logistic normal with added args for osa residuals
//do_mult = 1: do multiplicative transformation rather than additive
template<class Type>
Type dlogisticnormal_osa(vector<Type> obs, vector<Type> mu,  matrix<Type> S, data_indicator<vector<Type>, Type> keep, int do_mult, int do_log)
{
  using namespace density;
  MVNORM_t<Type> mvnorm(S);

  vector<Type> x(obs.size()-1);
  if(do_mult == 1){
    x = log(obs.head(obs.size()-1));
    for(int i = 0; i < x.size(); i++) x(i) -= log(1-x.head(i+1).sum());
  }
  else x = log(obs.head(obs.size()-1)) - log(obs(obs.size()-1));
  Type nll = mvnorm(x-mu, keep.head(x.size()));
  if(sum(keep)>x.size()) nll += log(obs).sum(); //jacobian, do it only when osa residuals not being calculated?

  if(do_log == 1) return -nll;
  else return exp(-nll);
}


template<class Type>
Type get_acomp_ll_osa(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref, data_indicator<vector<Type>, Type> keep)
{
  Type ll = 0.0;
  if(age_comp_model == 1) //multinomial
  {
    vector<Type> p = paa_pred + 1.0e-15;
    vector<Type> x = Neff * paa_obs;
    ll = dmultinom_osa(x, p, keep, 1);
  }
  if(age_comp_model == 2) //dirichlet-multinomial
  {
    vector<Type> x = Neff * paa_obs;
    vector<Type> p = paa_pred + 1.0e-15;
    vector<Type> alphas = p * exp(age_comp_pars(0));
    ll = ddirmultinom_osa(x, p, alphas, keep, 1);
  }
  if(age_comp_model == 3) //dirichlet
  {
    Type obs = 0.0, pred = 0.0, obs_2 = 0.0, pred_2 = 0.0;
    int npos = 0;
    for(int a = 0; a < paa_obs.size(); a++)
    {
      if(paa_obs(a) > Type(1.0e-15)) npos++;
      pred += paa_pred(a);
      obs += paa_obs(a);
      if(paa_obs(a) > Type(1.0e-15))

    for(int a = aref-1; a < n_ages; a++)
    {
      obs_2 += paa_obs(a);
      pred_2 += paa_pred(a);
    }
    ll = lgamma(exp(age_comp_pars(0)));
    for(int a = 0; a < aref-1; a++)
    {
      pred += paa_pred(a);
      obs += paa_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        ll +=  t_keep(a) * (-lgamma(exp(age_comp_pars(0)) * pred) + (exp(age_comp_pars(0)) * pred - 1.0) * log(obs));
        pred = 0.0;
        obs = 0.0;
      }
      //else pooling with next age
    }
    //add in the last age class(es).
    ll += t_keep(aref-1) * (-lgamma(exp(age_comp_pars(0)) * pred_2) + (exp(age_comp_pars(0)) * pred_2 - 1.0) * log(obs_2));
  }
  if(age_comp_model == 4) //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  {
    vector<Type> X(n_ages), p0(n_ages);
    Type mu = 0.0, sd = 0.0, pos_obs = 0.0, pos_pred = 0.0, pos_obs_l = 0.0, pos_pred_l = 0.0, pos_obs_sum = 0.0;
    Type pos_pred_sum = 0.0, y = 0.0;
    X = log(paa_pred + Type(1.0e-15)) - log(1.0 - paa_pred + Type(1.0e-15));
    p0 = 1.0/(1.0 + exp(exp(age_comp_pars(1))*(X - age_comp_pars(0)))); //prob of zero declines with proportion caught
    sd = exp(age_comp_pars(2));
    int last_pos = 0;
    pos_obs_sum = sum(paa_obs);
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
    {
      pos_pred_sum += paa_pred(a);
      last_pos = a;
    }
    //logistic applies only to proportions of non-zero observations
    pos_obs_l = paa_obs(last_pos)/pos_obs_sum;
    pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
    for(int a = 0; a < n_ages; a++)
    {
      if(paa_obs(a) < Type(1.0e-15)) ll += (log(p0(a) + Type(1.0e-15)));
      // if(paa_obs(a) < Type(1.0e-15)) ll += t_keep(a) * (log(p0(a) + Type(1.0e-15)));
      else
      {
        ll += (log(1.0 - p0(a) + Type(1.0e-15)));
        // ll += t_keep(a) * (log(1.0 - p0(a) + Type(1.0e-15)));
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += t_keep(a) * (-0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs));
        }
      }
    }
    ll -= t_keep(last_pos) * log(pos_obs_l); //add in the last observed age class(es).
  }
  if(age_comp_model == 5) //logistic normal. Pool zero observations with adjacent age classes.
  {
    Type mu = 0.0, sd = 0.0, obs = 0.0, pred = 0.0, obs_2 = 0.0, pred_2 = 0.0, y = 0.0;
    for(int a = aref-1; a < n_ages; a++)
    {
      obs_2 += paa_obs(a);
      pred_2 += paa_pred(a);
    }
    for(int a = 0; a < aref-1; a++)
    {
      pred += paa_pred(a);
      obs += paa_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        sd = exp(age_comp_pars(0)-0.5*log(Neff));
        y = log(obs) - log(obs_2);
        mu = log(pred + Type(1.0e-15)) - log(pred_2 + Type(1.0e-15));
        ll += t_keep(a) * (-0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(obs));
        pred = 0.0;
        obs = 0.0;
      }
      //else pooling with next age
    }
    for(int a = aref-1; a < n_ages; a++)
    {
      ll -= t_keep(a) * log(obs_2); //add in the last age class(es).
    }
  }
  if(age_comp_model == 6) //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
  {
    vector<Type> p0(n_ages);
    Type n_e = 0.0, mu = 0.0, sd = 0.0, pos_obs = 0.0, pos_pred = 0.0, pos_obs_l = 0.0, pos_pred_l = 0.0, pos_obs_sum = 0.0;
    Type pos_pred_sum = 0.0, y = 0.0;
    n_e = exp(age_comp_pars(0));
    p0 = exp(n_e * log(1.0-paa_pred + Type(1.0e-15))); //prob of zero declines with proportion caught
    sd = exp(age_comp_pars(1));
    int last_pos = 0;
    pos_obs_sum = sum(paa_obs);
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
    {
      pos_pred_sum += paa_pred(a);
      last_pos = a;
    }
    //logistic applies only to proportions of non-zero observations
    pos_obs_l = paa_obs(last_pos)/pos_obs_sum;
    pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
    for(int a = 0; a < n_ages; a++)
    {
      if(paa_obs(a) < Type(1.0e-15)) ll += (log(p0(a) + Type(1.0e-15)));
      // if(paa_obs(a) < Type(1.0e-15)) ll += t_keep(a) * (log(p0(a) + Type(1.0e-15)));
      else
      {
        ll += (log(1.0 - p0(a) + Type(1.0e-15)));
        // ll += t_keep(a) * (log(1.0 - p0(a) + Type(1.0e-15)));
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += t_keep(a) * (-0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs));
        }
      }
    }
    ll -= t_keep(last_pos) * log(pos_obs_l); //add in the last observed age class(es).
  }
  if(age_comp_model == 7) //logistic normal treating 0 observations as missing. One parameter.
  {
    Type mu = 0.0, pos_obs = 0.0, pos_pred = 0.0, pos_obs_l = 0.0, pos_pred_l = 0.0, pos_obs_sum = 0.0;
    Type pos_pred_sum = 0.0, y = 0.0;
    Type sd = exp(age_comp_pars(0));
    int last_pos = 0;
    pos_obs_sum = sum(paa_obs);
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
    {
      pos_pred_sum += paa_pred(a);
      last_pos = a;
    }
    //logistic applies only to proportions of non-zero observations
    pos_obs_l = paa_obs(last_pos)/pos_obs_sum;
    pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
    for(int a = 0; a < n_ages; a++)
    {
      if(paa_obs(a) > Type(1.0e-15))
      {
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += t_keep(a) * (-0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs));
        }
      }
    }
    ll -= log(pos_obs_l); //add in the last observed age class(es).
  }
  return ll;
}
