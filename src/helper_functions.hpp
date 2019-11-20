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
    if(obs(a)>0) ll += obs(a) * log(pred(a));
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
    if(obs(a) > 0) ll += t_keep(a) * (-lgamma(Neff*obs(a) + 1.0) + Neff*obs(a) * log(pred(a)));
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
  for(int i = 0; i < n; i++) ll +=  -lgamma(phi * p(i)) + (phi * p(i) - 1.0) * log(obs(i));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> p,  Type phi, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + phi * p(a)) - lgamma(phi * p(a));
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
        ll +=  -lgamma(exp(age_comp_pars(0)) * pred) + (exp(age_comp_pars(0)) * pred - 1.0) * log(obs);
        pred = 0.0;
        obs = 0.0;
      }
      //else pooling with next age
    }
    //add in the last age class(es).
    ll += -lgamma(exp(age_comp_pars(0)) * pred_2) + (exp(age_comp_pars(0)) * pred_2 - 1.0) * log(obs_2);
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
vector<Type> sim_acomp(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
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
    Type sd = exp(age_comp_pars(0));

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
matrix<Type> get_SPR_res(matrix<Type> MAA, matrix<Type> FAA, int which_F_age, array<Type> waa, int waa_pointer_ssb, int waa_pointer_tot_catch,
  matrix<Type> mature, Type percentSPR, vector<Type> predR, vector<Type> fracyr_SSB, vector<Type> log_SPR0)
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
  log_FXSPR_iter.fill(log(0.2)); //starting value
  res.col(5).fill(log(0.2));
  vector<Type> log_YPR_FXSPR(ny), sel(na), waacatch(na), waassb(na), mat(na), M(na);
  for(int y = 0; y < ny; y++)
  {
    M = MAA.row(y);
    waassb = ssbWAA.row(y);
    waacatch = catchWAA.row(y);
    mat = mature.row(y);
    sel = FAA.row(y)/FAA(y,which_F_age-1);
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

template <class Type>
Type get_F_from_Catch(Type Catch, vector<Type> NAA, vector<Type> M, vector<Type> sel, vector<Type> waacatch)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(0.2)); //starting value
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
Type get_FXSPR(vector<Type> M, vector<Type> sel, vector<Type> waacatch, vector<Type> waassb,
  vector<Type> mat, Type percentSPR, Type fracyr_SSB, Type log_SPR0)
{
  int n = 10;
  vector<Type> log_FXSPR_i(1);
  vector<Type> log_FXSPR_iter(n);
  log_FXSPR_iter.fill(log(0.2)); //starting value
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
