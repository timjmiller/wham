#ifndef WHAM_REFIT_LEGACY_HPP
#define WHAM_REFIT_LEGACY_HPP

inline vector<int> normalize_legacy_recruitment(vector<int>& recruit_model,
  int n_stocks, int n_regions, int n_seasons)
{
  vector<int> legacy_steepness(n_stocks);
  legacy_steepness.setZero();
  for(int stock = 0; stock < n_stocks; stock++) {
    if(recruit_model(stock) == 5 || recruit_model(stock) == 6) {
      if(n_stocks != 1 || n_regions != 1 || n_seasons != 1)
        error("Legacy steepness requires one stock, region and season.");
      legacy_steepness(stock) = recruit_model(stock) - 2;
      recruit_model(stock) -= 2;
    }
  }
  return legacy_steepness;
}

template<class Type>
void apply_legacy_steepness(const vector<int>& legacy_steepness,
  const matrix<Type>& mean_rec_pars, array<Type>& log_M,
  array<Type>& mature, array<Type>& waa_ssb,
  const matrix<Type>& fracyr_SSB, int n_years_model, int n_years_proj,
  matrix<Type>& log_SR_a, matrix<Type>& log_SR_b)
{
  int n_ages = log_M.dim(3);
  for(int stock = 0; stock < legacy_steepness.size(); stock++) if(legacy_steepness(stock)) {
    if(n_years_proj > 0) error("Legacy steepness does not support stored projections.");
    Type steepness = legacy_steepness(stock) == 3 ?
      Type(0.2) + Type(0.8) / (Type(1) + exp(-mean_rec_pars(stock,0))) :
      Type(0.2) + exp(mean_rec_pars(stock,0));
    for(int year = 0; year < n_years_model; year++) {
      Type survival = 1;
      Type spr = 0;
      for(int age = 0; age < n_ages; age++) {
        Type mortality = exp(log_M(stock,0,year,age));
        Type abundance = survival;
        if(age == n_ages - 1) abundance /= Type(1) - exp(-mortality);
        spr += abundance * mature(stock,year,age) * waa_ssb(stock,year,age) *
          exp(-fracyr_SSB(year,stock) * mortality);
        survival *= exp(-mortality);
      }
      Type log_alpha;
      Type log_beta;
      if(legacy_steepness(stock) == 3) {
        log_alpha = log(Type(4) * steepness / ((Type(1) - steepness) * spr));
        log_beta = log((Type(5) * steepness - Type(1)) /
          ((Type(1) - steepness) * spr)) - mean_rec_pars(stock,1);
      } else {
        Type productivity = Type(1.25) * log(Type(5) * steepness);
        log_alpha = productivity - log(spr);
        log_beta = log(productivity / spr) - mean_rec_pars(stock,1);
      }
      log_SR_a(year,stock) += log_alpha - mean_rec_pars(stock,0);
      log_SR_b(year,stock) += log_beta - mean_rec_pars(stock,1);
    }
  }
}

#endif