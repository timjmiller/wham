template <class Type>
Type square(Type x){return x*x;}
#define see(object) std::cout << #object ":\n" << object << "\n";

// transformation to ensure correlation parameters are between -1 and 1
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
matrix<Type> extract_matrix_array3(array<Type> a, int index) //matrix has to be the last two dims of the 3d array
{
  int nr = a.dim(1);
  int nc = a.dim(2);
  matrix<Type> mat(nr, nc);
  for(int i = 0; i < nr; i++) for(int j = 0; j < nc; j++) mat(i,j) = a(index,i,j);
  return mat;
}


template<class Type>
Type mydmultinom(vector<Type> obs, vector<Type> pred, int do_log)
{
  //multinomial
  int dim = obs.size();
  Type N = obs.sum();
  Type ll = lgamma(N + 1.0);
  for(int a = 0; a < dim; a++)
  {
    ll += -lgamma(obs(a) + 1.0);
    if(obs(a)>0) ll += obs(a) * log(pred(a) + Type(1.0e-15));
  }
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
Type mydmultinom_osa(vector<Type> obs, vector<Type> pred, Type Neff, int do_log, vector<Type> t_keep)
{
  //multinomial
  int dim = obs.size();
  Type N = Neff * obs.sum();
  Type ll = lgamma(N + 1.0);
  for(int a = 0; a < dim; a++)
  {
    if(obs(a) <= 0) ll += t_keep(a) * -lgamma(Neff*obs(a) + 1.0);
    if(obs(a) > 0) ll += t_keep(a) * (-lgamma(Neff*obs(a) + 1.0) + Neff*obs(a) * log(pred(a) + Type(1.0e-15)));
    // if(obs(a)>0) ll += t_keep(a)* obs(a) * log(pred(a));
  }
  if(do_log == 1) return ll;
  else return exp(ll);
}
/*
template<class Type>
Type mydmultinom(vector<Type> obs, vector<Type> pred, int do_log)
{
  //multinomial
  int dim = obs.size();
  Type N = obs.sum();
  Type ll = lgamma(N + 1.0);
  Type tot_pred = 0.0;
  for(int a = 0; a < dim; a++)
  {
    ll += -lgamma(obs(a) + 1.0);
    if(obs(a)>0)
    {
      ll += obs(a) * log(pred(a));
      tot_pred += pred(a);
    }
  }
  see(tot_pred);
  ll -= N * log(tot_pred);
  if(do_log == 1) return ll;
  else return exp(ll);
}
*/
template<class Type>
vector<Type> rmultinom(Type N, vector<Type> p)
{
  //multinomial
  int dim = p.size();
  vector<Type> x(dim);
  int Nint = CppAD::Integer(N);
  x.setZero();
  for(int i = 0; i < Nint; i++)
  {
    Type y = runif(0.0,1.0);
    for(int a = 0; a < dim; a++) if(y < p.head(a+1).sum())
    {
      x(a) += 1.0;
      break;
    }
  }
  return x;
}
template<class Type>
vector<Type> rdirichlet(vector<Type> p, Type phi)
{
  vector<Type> alpha = p * phi;
  vector<Type> obs = rgamma(alpha,Type(1.0));
  obs = obs/obs.sum();
  return obs;
}
template <class Type>
Type ddirichlet(vector<Type> obs, vector<Type>p, Type phi, int do_log)
{
  int n = obs.size();
  Type ll = lgamma(phi);
  for(int i = 0; i < n; i++) ll +=  -lgamma(phi * (p(i) + Type(1.0e-15))) + (phi * (p(i) + Type(1.0e-15)) - 1.0) * log(obs(i));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> p,  Type phi, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + phi * (p(a) + Type(1.0e-15))) - lgamma(phi * (p(a) + Type(1.0e-15)));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type> // modified for osa residuals
Type ddirmultinom_osa(vector<Type> obs, vector<Type> p,  Type phi, int do_log, vector<Type> t_keep)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += t_keep(a) * (-lgamma(obs(a) + 1.0) + lgamma(obs(a) + phi * p(a)) - lgamma(phi * p(a)));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
vector<Type> rdirmultinom(Type N, vector<Type> p, Type phi) //dirichlet generated from iid gammas
{
  int Nint = CppAD::Integer(N);
  int dim = p.size();
  vector<Type> obs(dim);
  obs.setZero();
  for(int i = 0; i < Nint; i++)
  {
    vector<Type> dp = rdirichlet(p, phi);
    obs = obs + rmultinom(Type(1),dp);
  }
  return(obs);
}

template<class Type>
Type get_acomp_ll(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  Type ll = 0.0;
  if(age_comp_model == 1) //multinomial
  {
    vector<Type> temp_n = Neff * paa_obs;
    ll = mydmultinom(temp_n, paa_pred, 1);
  }
  if(age_comp_model == 2) //dirichlet-multinomial
  {
    vector<Type> temp_n = Neff * paa_obs;
    ll = ddirmultinom(temp_n, paa_pred, exp(age_comp_pars(0)),1);
  }
  if(age_comp_model == 3) //dirichlet
  {
    Type obs = 0.0, pred = 0.0, obs_2 = 0.0, pred_2 = 0.0;
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
        ll +=  -lgamma(exp(age_comp_pars(0)) * (pred + Type(1.0e-15))) + (exp(age_comp_pars(0)) * (pred + Type(1.0e-15)) - 1.0) * log(obs);
        pred = 0.0;
        obs = 0.0;
      }
      //else pooling with next age
    }
    //add in the last age class(es).
    ll += -lgamma(exp(age_comp_pars(0)) * (pred_2 + Type(1.0e-15))) + (exp(age_comp_pars(0)) * (pred_2 + Type(1.0e-15)) - 1.0) * log(obs_2);
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
      if(paa_obs(a) < Type(1.0e-15)) ll += log(p0(a) + Type(1.0e-15));
      else
      {
        ll += log(1.0 - p0(a) + Type(1.0e-15));
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += -0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
        }
      }
    }
    ll -= log(pos_obs_l); //add in the last observed age class(es).
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
        ll += -0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(obs);
        pred = 0.0;
        obs = 0.0;
      }
      //else pooling with next age
    }
    ll -= log(obs_2); //add in the last age class(es).
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
      if(paa_obs(a) < Type(1.0e-15)) ll += log(p0(a) + Type(1.0e-15));
      else
      {
        ll += log(1.0 - p0(a) + Type(1.0e-15));
        if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
        {
          pos_pred = paa_pred(a)/pos_pred_sum;
          pos_obs = paa_obs(a)/pos_obs_sum;
          y = log(pos_obs) - log(pos_obs_l);
          mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
          ll += -0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
        }
      }
    }
    ll -= log(pos_obs_l); //add in the last observed age class(es).
  }
  if(age_comp_model == 7) //logistic normal treating 0 observations as missing. One parameter.
  {
    Type mu = 0.0, pos_obs = 0.0, pos_pred = 0.0, pos_obs_l = 0.0, pos_pred_l = 0.0, pos_obs_sum = 0.0;
    Type pos_pred_sum = 0.0, y = 0.0;
    Type sd = exp(age_comp_pars(0)-0.5*log(Neff));
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
          ll += -0.5 * (log(2.0 * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
        }
      }
    }
    ll -= log(pos_obs_l); //add in the last observed age class(es).
  }
  return ll;
}

template<class Type>
Type get_acomp_ll_osa(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref, vector<Type> t_keep)
{
  Type ll = 0.0;
  if(age_comp_model == 1) //multinomial
  {
    // vector<Type> temp_n = Neff * paa_obs;
    ll = mydmultinom_osa(paa_obs, paa_pred, Neff, 1, t_keep);
    // vector<Type> temp_n = Neff * paa_obs;
    // ll = mydmultinom_osa(temp_n, paa_pred, 1, t_keep);
  }
  if(age_comp_model == 2) //dirichlet-multinomial
  {
    vector<Type> temp_n = Neff * paa_obs;
    ll = ddirmultinom_osa(temp_n, paa_pred, exp(age_comp_pars(0)),1, t_keep);
  }
  if(age_comp_model == 3) //dirichlet
  {
    Type obs = 0.0, pred = 0.0, obs_2 = 0.0, pred_2 = 0.0;
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

template<class Type>
vector<Type> sim_acomp(Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  int n_ages = paa_obs.size();
  vector<Type> obs(n_ages);
  obs.setZero();
  if(age_comp_model == 1)
  {
    //int N = CppAD::Integer(Neff);
    obs = rmultinom(Neff, paa_pred);
    obs = obs/obs.sum();// proportions
  }
  if(age_comp_model == 2) //dirichlet-multinomial. dirichlet generated from iid gammas and multinomial from uniform
  {
    //int N = CppAD::Integer(Neff);
    obs = rdirmultinom(Neff,paa_pred,exp(age_comp_pars(0)));
    obs = obs/obs.sum();// proportions
  }
  if(age_comp_model == 3) //dirichlet generated from iid gammas
  {
    Type obs_2 = 0.0;
    vector<Type> best_obs = rdirichlet(paa_pred, exp(age_comp_pars(0)));
    obs_2 = best_obs.tail(n_ages-aref+1).sum(); // .tail last n_ages-aref+1 components
    for(int a = aref-1; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15)) obs(a) = obs_2;
    obs_2 = 0.0;
    for(int a = 0; a < aref-1; a++)
    {
      obs_2 += best_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        obs(a) = obs_2;
        obs_2 = 0.0;
      }
      else obs(a) = 0.0;
      //else pooling with next age
    }
  }
  if(age_comp_model == 4) //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
  {
    vector<Type> X = log(paa_pred + Type(1.0e-15)) - log(1.0 - paa_pred + Type(1.0e-15));
    vector<Type> p0 = 1.0/(1.0 + exp(exp(age_comp_pars(1))*(X - age_comp_pars(0)))); //prob of zero declines with proportion caught
    Type sd = exp(age_comp_pars(2));
    for(int a = 0; a < n_ages; a++) obs(a) = rbinom(Type(1.0), Type(1.0) - p0(a)); // generate instances of positive observations
    int n_pos = 0;
    for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5) n_pos++;
    if(n_pos>0)
    {
      vector<Type> pos_pred(n_pos);
      int k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        pos_pred(k) = paa_pred(a);
        k++;
      }
      vector<Type> pos_obs(n_pos);
      pos_obs.setZero();
      for(int a = 0; a < n_pos-1; a++)
      {
        pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
      }
      pos_obs = pos_obs/(1.0 + pos_obs.sum());
      pos_obs(n_pos-1) = 1.0 - pos_obs.sum();
      k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        obs(a) = pos_obs(k);
        k++;
      }
    }
  }
  if(age_comp_model == 5) //logistic normal. Pool zero observations with adjacent age classes.
  {
    vector<Type> best_obs(n_ages);
    best_obs.setZero();
    Type sd = exp(age_comp_pars(0)-0.5*log(Neff));
    for(int a = 0; a < n_ages-1; a++) best_obs(a) = exp(rnorm(log(paa_pred(a)) - log(paa_pred(n_ages-1)), sd));
    best_obs = best_obs/(1.0 + best_obs.sum());
    best_obs(n_ages-1) = 1.0 - best_obs.sum();

    Type obs_2 = best_obs.tail(n_ages-aref+1).sum(); // .tail last n_ages-aref+1 components
    for(int a = aref-1; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15)) obs(a) = obs_2;
    obs_2 = 0.0;
    for(int a = 0; a < aref-1; a++)
    {
      obs_2 += best_obs(a);
      if(paa_obs(a) > Type(1.0e-15))
      {
        obs(a) = obs_2;
        obs_2 = 0.0;
      }
      else obs(a) = 0.0;
      //else pooling with next age
    }
  }
  if(age_comp_model == 6) //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
  {
    Type n_e = exp(age_comp_pars(0));
    vector<Type> p0 = exp(n_e * log(1.0-paa_pred + Type(1.0e-15))); //prob of zero declines with proportion caught
    Type sd = exp(age_comp_pars(1));

    for(int a = 0; a < n_ages; a++) obs(a) = rbinom(Type(1.0), Type(1.0) - p0(a)); // generate instances of positive observations
    int n_pos = 0;
    for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5) n_pos++;
    if(n_pos>0)
    {
      vector<Type> pos_pred(n_pos);
      int k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        pos_pred(k) = paa_pred(a);
        k++;
      }
      vector<Type> pos_obs(n_pos);
      pos_obs.setZero();
      for(int a = 0; a < n_pos-1; a++)
      {
        pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
      }
      pos_obs = pos_obs/(1.0 + pos_obs.sum());
      pos_obs(n_pos-1) = 1.0 - pos_obs.sum();
      k = 0;
      for(int a = 0; a < n_ages; a++) if(obs(a) > 0.5)
      {
        obs(a) = pos_obs(k);
        k++;
      }
    }
  }
  if(age_comp_model == 7) //logistic normal treating 0 observations as missing. One parameter.
  {
    Type sd = exp(age_comp_pars(0)-0.5*log(Neff));

    int n_pos = 0;
    for(int a = 0; a < n_ages; a++) if(paa_obs(a) > 1.0e-15) n_pos++;
    if(n_pos>0)
    {
      vector<Type> pos_pred(n_pos);
      int k = 0;
      for(int a = 0; a < n_ages; a++) if(paa_obs(a) > 1.0e-15)
      {
        pos_pred(k) = paa_pred(a);
        k++;
      }
      vector<Type> pos_obs(n_pos);
      pos_obs.setZero();
      for(int a = 0; a < n_pos-1; a++)
      {
        pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
      }
      pos_obs = pos_obs/(1.0 + pos_obs.sum());
      pos_obs(n_pos-1) = 1.0 - pos_obs.sum();
      k = 0;
      for(int a = 0; a < n_ages; a++) if(paa_obs(a) > 1.0e-15)
      {
        obs(a) = pos_obs(k);
        k++;
      }
    }
  }
  return obs;
}

template <class Type>
vector<matrix<Type> > get_selectivity(int n_years, int n_ages, int n_selblocks, vector<matrix<Type> > selpars, vector<int> selblock_models)
{
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

template <class Type>
Type get_SPR(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type SPR = 0.0, ntemp = 1.0;
  vector<Type> F(n_ages), Z(n_ages);

  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    SPR += ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0-exp(-Z(n_ages-1));
  SPR += ntemp * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB*Z(n_ages-1));

  return SPR;
}

template <class Type>
vector<Type> get_SPRAA(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type ntemp = 1.0;
  vector<Type> F(n_ages), Z(n_ages), SPR(n_ages);

  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    SPR(age) = ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0-exp(-Z(n_ages-1));
  SPR(n_ages-1) = ntemp * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB*Z(n_ages-1));

  return SPR;
}

template <class Type>
Type get_YPR(Type log_F, vector<Type> M, vector<Type> sel, vector<Type> waacatch)
{
  int n_ages = M.size();
  Type YPR = 0.0, ntemp = 1.0;
  vector<Type> F(n_ages), Z(n_ages);

  F = exp(log_F) * sel;
  Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    YPR += ntemp * F(age) * waacatch(age) * (1.0 - exp(-Z(age)))/Z(age);
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0 - exp(-Z(n_ages-1));
  YPR += ntemp * F(n_ages-1) * waacatch(n_ages-1) * (1.0 - exp(-Z(n_ages-1)))/Z(n_ages-1);

  return YPR;
}
template <class Type>
Type get_SPR_0(vector<Type> M, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  Type SPR_0 = Type(0.0);
  Type ntemp0 = Type(1.0);
  for (int a = 0; a < n_ages - 1; a++)
  {
    SPR_0 += ntemp0 * mat(a) * waassb(a) * exp(-fracyearSSB * M(a));
    ntemp0 *= exp(-M(a));
  }
  ntemp0 /= Type(1.0)-exp(-M(n_ages-1));
  SPR_0 += ntemp0 * mat(n_ages-1) * waassb(n_ages-1) * exp(-fracyearSSB * M(n_ages-1));
  return SPR_0;
}


/* calculate beverton-holt or ricker equilibrium yield */
template<class Type>
struct sr_yield {
  /* Data and parameter objects for yield calculation: */
  Type SR_a;
  Type SR_b;
  vector<Type> M;
  vector<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  vector<Type> waacatch;
  Type fracyearSSB;
  int sr_type;

  /* Constructor */
  sr_yield(Type SR_a_, Type SR_b_,
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  vector<Type> waacatch_,
  Type fracyearSSB_, int sr_type_) :
    SR_a(SR_a_), SR_b(SR_b_), M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), waacatch(waacatch_), fracyearSSB(fracyearSSB_), sr_type(sr_type_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    T YPR = 0, SPR = 0, ntemp = 1, R;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++)
    {
      YPR += ntemp * F(age) * T(waacatch(age)) * (1- exp(-Z(age)))/Z(age);
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    YPR += ntemp * F(n_ages-1) * T(waacatch(n_ages-1)) * (1 - exp(-Z(n_ages-1)))/Z(n_ages-1);
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    //Type SPR = get_SPR(x, M, sel, mat, waassb, fracyearSSB);
    if(sr_type == 0) R = (T(SR_a) - 1/SPR) / T(SR_b); //beverton-holt
    if(sr_type == 1) R = log(T(SR_a) * SPR)/(T(SR_b) * SPR); //ricker
    T Y = YPR * R;
    return Y;
  }
};

/* calculate SPR at F */
template<class Type>
struct spr_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  vector<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  Type fracyearSSB;

  /* Constructor */
  spr_F(
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  Type fracyearSSB_) :
    M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), fracyearSSB(fracyearSSB_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    T SPR = 0, ntemp = 1;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++)
    {
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    return SPR;
  }
};

template <class Type>
matrix<Type> get_SPR_res(matrix<Type> MAA, matrix<Type> FAA, vector<int> which_F_age, array<Type> waa, int waa_pointer_ssb, int waa_pointer_tot_catch,
  matrix<Type> mature, Type percentSPR, vector<Type> predR, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<Type> F_init)
{
  int n = 10;
  int ny = MAA.rows();
  int na = MAA.cols();
  matrix<Type> res(ny, 5+n);
  vector<Type> log_FXSPR(ny), log_FXSPR_i(1), log_SPR_FXSPR(ny);
  vector<Type> log_SSB_FXSPR(ny), log_Y_FXSPR(ny);
  matrix<Type> log_FXSPR_iter(ny,n);
  matrix<Type> catchWAA = extract_matrix_array3(waa, waa_pointer_tot_catch-1);
  matrix<Type> ssbWAA = extract_matrix_array3(waa, waa_pointer_ssb-1);
//  log_FXSPR_iter.fill(log(0.2)); //starting value
//  res.col(5).fill(log(0.2));
  log_FXSPR_iter.col(0) = log(F_init);
  res.col(5) = log_FXSPR_iter.col(0);
  vector<Type> log_YPR_FXSPR(ny), sel(na), waacatch(na), waassb(na), mat(na), M(na);
  for(int y = 0; y < ny; y++)
  {
    M = MAA.row(y);
    waassb = ssbWAA.row(y);
    waacatch = catchWAA.row(y);
    mat = mature.row(y);
    sel = FAA.row(y)/FAA(y,which_F_age(y)-1);
    spr_F<Type> sprF(M, sel, mat, waassb, fracyr_SSB(y));
    //log_SPR0(y) = log(get_SPR_0(M, mat, waassb, fracyr_SSB(y)));
    for (int i=0; i<n-1; i++)
    {
      log_FXSPR_i(0) = log_FXSPR_iter(y,i);
      vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
      //matrix<Type> hess_sr_yield = autodiff::hessian(sprF,log_FXSPR_i);
      log_FXSPR_iter(y,i+1) = log_FXSPR_iter(y,i) - (sprF(log_FXSPR_i) - 0.01*percentSPR*exp(log_SPR0(y)))/grad_spr_F(0);// /hess_sr_yield(0,0);
      res(y,5+i+1) = log_FXSPR_iter(y,i+1);
    }
    log_FXSPR(y) = log_FXSPR_iter(y,n-1);
    log_SPR_FXSPR(y) = log(get_SPR(log_FXSPR(y), M, sel, mat, waassb, fracyr_SSB(y)));
    log_YPR_FXSPR(y) = log(get_YPR(log_FXSPR(y), M, sel, waacatch));
    log_SSB_FXSPR(y) = log(predR(y)) + log_SPR_FXSPR(y);
    log_Y_FXSPR(y) = log(predR(y)) + log_YPR_FXSPR(y);
    res(y,0) = log_FXSPR(y);
    res(y,1) = log_SSB_FXSPR(y);
    res(y,2) = log_Y_FXSPR(y);
    //res(y,3) = log_SPR0(y);
    res(y,3) = log_SPR_FXSPR(y);
    res(y,4) = log_YPR_FXSPR(y);
  }
  return res;
}

/* calculate Catch at F */
template<class Type>
struct catch_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  vector<Type> sel;
  vector<Type> waacatch;
  vector<Type> NAA;

  /* Constructor */
  catch_F(
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> waacatch_,
  vector<Type> NAA_) :
    M(M_), sel(sel_), waacatch(waacatch_), NAA(NAA_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_ages = M.size();
    T Catch = 0;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++)
    {
      Catch += T(waacatch(age)) * T(NAA(age)) * F(age) *(1 - exp(-Z(age)))/Z(age);
    }
    return Catch;
  }
};

/* calculate Catch at F */
template<class Type>
struct log_catch_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  vector<Type> sel;
  vector<Type> waacatch;
  vector<Type> NAA;

  /* Constructor */
  log_catch_F(
  vector<Type> M_,
  vector<Type> sel_,
  vector<Type> waacatch_,
  vector<Type> NAA_) :
    M(M_), sel(sel_), waacatch(waacatch_), NAA(NAA_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_ages = M.size();
    T Catch = 0;
    vector<T> F(n_ages), Z(n_ages);

    F = exp(log_F(0)) * sel.template cast<T>();
    Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++)
    {
      Catch += T(waacatch(age)) * T(NAA(age)) * F(age) *(1 - exp(-Z(age)))/Z(age);
    }
    return log(Catch);
  }
};

template <class Type>
Type get_F_from_Catch(Type Catch, vector<Type> NAA, vector<Type> M, vector<Type> sel, vector<Type> waacatch, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  catch_F<Type> catchF(M, sel, waacatch, NAA);
  for (int i=0; i<n-1; i++)
  {
    log_F_i(0) = log_F_iter(i);
    vector<Type> grad_catch_F = autodiff::gradient(catchF,log_F_i);
    log_F_iter(i+1) = log_F_iter(i) - (catchF(log_F_i) - Catch)/grad_catch_F(0);
  }
  Type res = exp(log_F_iter(n-1));
  return res;
}

template <class Type>
Type get_F_from_log_Catch(Type Catch, vector<Type> NAA, vector<Type> M, vector<Type> sel, vector<Type> waacatch, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  log_catch_F<Type> logcatchF(M, sel, waacatch, NAA);
  for (int i=0; i<n-1; i++)
  {
    log_F_i(0) = log_F_iter(i);
    vector<Type> grad_log_catch_F = autodiff::gradient(logcatchF,log_F_i);
    log_F_iter(i+1) = log_F_iter(i) - (logcatchF(log_F_i) - log(Catch))/grad_log_catch_F(0);
  }
  Type res = exp(log_F_iter(n-1));
  return res;
}

template <class Type>
Type get_FXSPR(vector<Type> M, vector<Type> sel, vector<Type> waassb,
  vector<Type> mat, Type percentSPR, Type fracyr_SSB, Type log_SPR0, Type F_init)
{
  int n = 10;
  vector<Type> log_FXSPR_i(1);
  vector<Type> log_FXSPR_iter(n);
  log_FXSPR_iter.fill(log(F_init)); //starting value
  spr_F<Type> sprF(M, sel, mat, waassb, fracyr_SSB);
  for (int i=0; i<n-1; i++)
  {
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR*exp(log_SPR0))/grad_spr_F(0);// /hess_sr_yield(0,0);
  }
  Type res = exp(log_FXSPR_iter(n-1));
  return res;
}

template <class Type>
Type get_FMSY(Type log_a, Type log_b, vector<Type> M, vector<Type> sel, vector<Type> waacatch, vector<Type> waassb,
  vector<Type> mat, Type fracyr_SSB, Type log_SPR0, int recruit_model, Type F_init)
{    
  int n = 10;
  vector<Type> log_FMSY_i(1);
  vector<Type> log_FMSY_iter(n);
  log_FMSY_iter.fill(log(F_init)); //starting value
  Type a = exp(log_a);
  Type b = exp(log_b);
  int sr_type = 0; //recruit_model = 3, B-H
  if(recruit_model == 4) sr_type = 1; //recruit_model = 4, Ricker
  sr_yield<Type> sryield(a, b, M, sel, mat, waassb, waacatch, fracyr_SSB, sr_type);
  for (int i=0; i<n-1; i++)
  {
    log_FMSY_i(0) = log_FMSY_iter(i);
    vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
    matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
    log_FMSY_iter(i+1) = log_FMSY_iter(i) - grad_sr_yield(0)/hess_sr_yield(0,0);
  }
  Type FMSY = exp(log_FMSY_iter(n-1));
  return FMSY;
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
Type get_pred_recruit_y(int y, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSB, matrix<Type> NAA, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, vector<int> Ecov_where, vector<int> Ecov_how, matrix<Type> Ecov_lm){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: integer determining if Ecov is affecting recruitment
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: matrix that holds linear predictor for Ecov
   */
  //recruit_model == 1, random walk
  Type pred_recruit = NAA(y-1,0);
  if(recruit_model == 1) // random walk
  {
    //pred_NAA(y,0) = NAA(y-1,0);
  }
  else
  {
    if(recruit_model == 2) // random about mean
    {
      pred_recruit = exp(mean_rec_pars(0));
      int nE = Ecov_where.size();
      for(int i=0; i < nE; i++){
        if(Ecov_where(i) == 1) if(Ecov_how(i) == 1) pred_recruit *= exp(Ecov_lm(y,i));
      }
      //pred_NAA(y,0) = exp(mean_rec_pars(0));
      //if(Ecov_recruit > 0) if(Ecov_how(Ecov_recruit-1) == 1) pred_NAA(y,0) *= exp(Ecov_lm(y,Ecov_recruit-1));
    }
    else
    {
      if(recruit_model == 3) // BH stock recruit (if ecov effect, already modified SR_a and/or SR_b)
      {
        pred_recruit = exp(log_SR_a(y)) * SSB(y-1)/(1 + exp(log_SR_b(y))*SSB(y-1));
        //pred_NAA(y,0) = exp(log_SR_a(y)) * SSB(y-1)/(1 + exp(log_SR_b(y))*SSB(y-1));
      }
      else // recruit_model = 4, Ricker stock recruit (if ecov effect, already modified SR_a and/or SR_b)
      {
        pred_recruit = exp(log_SR_a(y)) * SSB(y-1) * exp(-exp(log_SR_b(y)) * SSB(y-1));
        //pred_NAA(y,0) = exp(log_SR_a(y)) * SSB(y-1) * exp(-exp(log_SR_b(y)) * SSB(y-1));
      }
    }
  }
  return(pred_recruit);
}

template <class Type>
vector<Type> get_pred_NAA_y(int y, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSB, matrix<Type> NAA, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, vector<int> Ecov_where, vector<int> Ecov_how, matrix<Type> Ecov_lm, matrix<Type> ZAA){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: integer vector determining if Ecov i affects recruitment (= 1)
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: matrix that holds linear predictor for Ecov
   * ZAA: matrix of total mortality rate by year and age
   */
  int n_ages = NAA.cols();
  vector<Type> pred_NAA(n_ages);
  
  // Expected recruitment
  pred_NAA(0) = get_pred_recruit_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
    log_SR_b, Ecov_where, Ecov_how, Ecov_lm);
    
  // calculate pred_NAA for ages after recruitment
  for(int a = 1; a < n_ages-1; a++) pred_NAA(a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
  pred_NAA(n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
  
  return(pred_NAA);
}

template <class Type>
vector<Type> get_waa_y(array<Type> waa, int y, int na, int pointer){
  vector<Type> waay(na);
  for(int a = 0; a < na; a++) waay(a) = waa(pointer-1, y, a);
  return(waay);
}

template <class Type>
Type get_SSB(matrix<Type> NAA, matrix<Type> ZAA, array<Type> waa, matrix<Type> mature, int y, int ssbpointer, vector<Type> fracyr_SSB){
  int na = mature.cols();
  Type SSB = 0;
  vector<Type> waay = get_waa_y(waa,y,na,ssbpointer);
  for(int a = 0; a < na; a++) SSB += NAA(y,a) * waay(a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
  return(SSB);
}
 
template <class Type>
matrix<Type> get_F_proj(int y, int n_fleets, vector<int> proj_F_opt, array<Type> FAA, matrix<Type> NAA, matrix<Type> MAA, matrix<Type> mature, 
  vector<Type> waacatch, vector<Type> waassb, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, vector<int> which_F_age, Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, Type F_init){
    /* 
     get F to project for next time step
              y:  year of projection (>n_years_model)
     proj_F_opt:  for each year, how to specify F for projection (1 to 5)
        FAA_tot:  FAA_tot matrix from main code.
            NAA:  NAA matrix from main code
            MAA:  MAA matrix from main code.
         mature:  maturity matrix from main code.
       waacatch:  weight at age in catch to use for year y (use function get_waa_y defined above)
         waassb:  weight at age for SSB to use in year y (use function get_waa_y defined above)
     fracyr_SSB:  vector of yearly fractions of the year when spawning occurs
       log_SPR0:  vector of yearly log(unfished SSB/R) 
  avg_years_ind:  integer vector of years to average F for projection
  n_years_model:  number of years before projection begins
    which_F_age:  define which age has max F
     percentSPR:  percentage (0-100) of unfished spawning potential to determine F_percentSPR
    proj_Fcatch: vector (n_years_proj) of user specified Fishing mortality rates to project
   percentFXSPR:  percentage (0-100) of F_percentSPR to use in catch, e.g. GOM cod uses F = 75% F_40%SPR
         F_init: initial value to use for FXSPR or FMSY newton method
    */
    //if(y > n_years_model-1){
    //for(int a = 0; a < n_ages; a++) waacatch(a) = waa(waa_pointer_totcatch-1, y, a);
    //for(int a = 0; a < n_ages; a++) waassb(a) = waa(waa_pointer_ssb-1, y, a);
  int n_toavg = avg_years_ind.size();
  int n_ages = waacatch.size();
  int proj_F_opt_y = proj_F_opt(y-n_years_model);

  //proj_F_opt == 1, last year F (default)
  matrix<Type> FAA_proj(n_fleets, n_ages);
  vector<Type> FAA_tot_proj(n_ages);
  FAA_tot_proj.setZero();
  if(proj_F_opt_y == 1){ // last year F (default)
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_proj(f,a) = FAA(n_years_model-1,f,a);
  }
  if(proj_F_opt_y>1)
  {
    //array<Type> FAA_toavg(n_fleets,n_toavg, n_ages);
    FAA_proj.setZero();
    //option 1: average F is by fleet and Ftot is sum of fleet averages
    //when there is more than 1 fleet, the sum of predicted catch across fleets will not generally equal the total catch using FAA_tot and waa_totcatch.
    for(int f = 0; f < n_fleets; f++) 
    {
      for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
        FAA_proj(f,a) += FAA(avg_years_ind(i),f,a);
      }
      FAA_proj.row(f) /= Type(n_toavg);
    }
    //get selectivity using average over avg.yrs
    matrix<Type> sel_proj(n_fleets,n_ages);
    vector<Type> sel_tot_proj(n_ages);
    vector<Type> FAA_tot_proj = FAA_proj.colwise().sum();
    for(int a = 0; a < n_ages; a++) sel_tot_proj(a) = FAA_tot_proj(a)/FAA_tot_proj(which_F_age(y)-1);
    //selectivity at age and fleet as a proportion of full total F (sum across fleets and ages = 1).
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) sel_proj(f,a) = FAA_proj(f,a)/FAA_tot_proj(which_F_age(y)-1);
    vector<Type> M = MAA.row(y);
    vector<Type> NAA_y = NAA.row(y);
    vector<Type> mat = mature.row(y);
    
    //proj_F_opt_y == 2, average F
    if(proj_F_opt_y == 2){
      //already defined
      //for(int f = 0; f < n_fleets; f++) FAA_proj = FAA_toavg.colwise().mean();
    }

    if(proj_F_opt_y == 3){ // F at X% SPR
      Type fracyr_SSB_y = fracyr_SSB(y);
      Type log_SPR0_y = log_SPR0(y);
      //FAA_tot.row(y) = get_FXSPR(M, sel_proj, waacatch, waassb, mat, percentSPR, fracyr_SSB_y, log_SPR0_y) * sel_proj;
      FAA_tot_proj = get_FXSPR(M, sel_tot_proj, waassb, mat, percentSPR, fracyr_SSB_y, log_SPR0_y, F_init) * sel_tot_proj;
      FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj * 0.01*percentFXSPR;
    }
    if(proj_F_opt_y == 4){ // user-specified F
      if(proj_Fcatch(y-n_years_model) < 1e-10){ // if F = 0, sel_proj is NaN
        FAA_proj.setZero();
      } else {
        //FAA_tot.row(y) = Type(proj_Fcatch(y-n_years_model)) * sel_proj;
        FAA_tot_proj = Type(proj_Fcatch(y-n_years_model)) * sel_tot_proj;
        FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
      }
    }
    if(proj_F_opt_y == 5){ // calculate F from user-specified catch
      Type thecatch = proj_Fcatch(y-n_years_model);
      if(thecatch < 1e-10){ // if catch = 0, F = 0 and sel_proj is NaN
        FAA_proj.setZero();
      } else {
        //FAA_tot.row(y) = get_F_from_Catch(thecatch, NAA_y, M, sel_proj, waacatch) * sel_proj;
        //FAA_tot_proj = get_F_from_Catch(thecatch, NAA_y, M, sel_tot_proj, waacatch) * sel_tot_proj;
        FAA_tot_proj = get_F_from_log_Catch(thecatch, NAA_y, M, sel_tot_proj, waacatch, F_init) * sel_tot_proj;
        FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
      }
    }
  }
  return(FAA_proj);
}


template <class Type>
matrix<Type> sim_pop(array<Type> NAA_devs, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSBin, matrix<Type> NAAin, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, vector<int> Ecov_where, vector<int> Ecov_how, matrix<Type> Ecov_lm, int n_NAA_sigma, 
  int do_proj, vector<int> proj_F_opt, array<Type> FAA, matrix<Type> FAA_tot, matrix<Type> MAA, matrix<Type> mature, array<Type> waa, 
  int waa_pointer_totcatch, int waa_pointer_ssb, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, int n_fleets, vector<int> which_F_age, Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, vector<Type> F_proj_init){

  // Population model (get NAA, numbers-at-age, for all years)
  int ny = log_SR_a.size();
  int n_ages = NAAin.cols();
  vector<Type> SSB(ny);
  SSB(0) = SSBin(0);
  matrix<Type> NAA(ny, n_ages), pred_NAA(ny, n_ages);
  matrix<Type> ZAA = MAA + FAA_tot;
  matrix<Type> log_NAA(ny-1, n_ages);
  log_NAA.setZero();
  NAA.row(0) = pred_NAA.row(0) = NAAin.row(0);
  for(int y = 1; y < ny; y++)
  {
    //use NAA.row(y-1)
    pred_NAA.row(y) = get_pred_NAA_y(y, recruit_model, mean_rec_pars, SSB, NAA, log_SR_a, 
      log_SR_b, Ecov_where, Ecov_how, Ecov_lm, ZAA);
    
    // calculate NAA
    if(n_NAA_sigma > 1){
      // all NAA are estimated (random effects)
      for(int a = 0; a < n_ages; a++) 
      {
        log_NAA(y-1,a) = log(pred_NAA(y,a)) + NAA_devs(y-1,a);
        NAA(y,a) = exp(log_NAA(y-1,a));
      }
    } else {
      // only recruitment estimated (either fixed or random effects)
      log_NAA(y-1,0) = log(pred_NAA(y,0)) + NAA_devs(y-1,0);
      NAA(y,0) = exp(log_NAA(y-1,0));
      for(int a = 1; a < n_ages; a++) NAA(y,a) = pred_NAA(y,a); // survival is deterministic     
    }

    // calculate F and Z in projection years, here bc need NAA(y) if using F from catch
    if(do_proj == 1){ // only need FAA_tot for projections, use average FAA_tot over avg.yrs
      // get selectivity using average over avg.yrs
      if(y > n_years_model-1){
        vector<Type> waacatch = get_waa_y(waa, y, n_ages, waa_pointer_totcatch);
        vector<Type> waassb = get_waa_y(waa, y, n_ages, waa_pointer_ssb);
        matrix<Type> FAA_proj = get_F_proj(y, n_fleets, proj_F_opt, FAA, NAA, MAA, mature, waacatch, waassb, fracyr_SSB, log_SPR0, avg_years_ind, n_years_model,
         which_F_age, percentSPR, proj_Fcatch, percentFXSPR, F_proj_init(y-n_years_model));
        for(int f = 0; f < n_fleets; f++) for(int a = 0; a< n_ages; a++) FAA(y,f,a) = FAA_proj(f,a);
        FAA_tot.row(y) = FAA_proj.colwise().sum();
        ZAA.row(y) = FAA_tot.row(y) + MAA.row(y);
      }
    } // end proj F
    SSB(y) = get_SSB(NAA,ZAA,waa, mature,y, waa_pointer_ssb, fracyr_SSB);
  } // end pop model loop
  
  matrix<Type> out(ny, 2 * n_ages + n_fleets * n_ages + 1);
  for(int i = 0; i < n_ages; i++) out.col(i) = NAA.col(i);
  for(int i = n_ages; i < 2 * n_ages; i++) out.col(i) = pred_NAA.col(i-n_ages);
  for(int y = 0; y < ny; y++) for(int i = 2 ; i < (2 + n_fleets); i++) for(int a = 0; a < n_ages; a++) out(y,i*n_ages + a) = FAA(y,i-2,a);
  out.col(out.cols()-1) = SSB;
  return(out);
}
