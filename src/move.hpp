template <class Type>
additive_ln_transform(vector<Type> x, int region, vector<int> do_move){
  //use additive transformation (e.g., logistic-normal model)
  //ensures that probabilities of moving and staying add to 1
  int D = x.size()+1;
  vector<Type> y(D);
  y.setZero();
  int j = 0;
  for(int i = 0; i < D; i++) {
    if(i != region) {
      if(do_move(i)==1) y(i) = exp(x(j)); //else prob of moving will be 0.
      j++;
    } else { //prob of staying will be 1- prob of moving
      y(i) = 1.0;
    }
  }
  y /= sum(y);
  return(y);
}


nll_mu_prior_re

simulate_mu_prior_re
  
nll_mu_re
  

simulate_mu_re
  
array<Type> get_mu(array<Type> trans_mu, array<Type> mu_re, array<int> can_move, int mig_type)
{
  if(can_move.sum()>0) //migration is happening
  {
    if(mig_type == 0) //migration is instantaneous after survival and mortality, so P is easy.
    {
      matrix<Type> move(n_regions,n_regions);
      move.setZero();
      //additive logistic-normal transform for probability of movement out of regions only; prob staying is 1 - sum(prob move)
      //for(int i = 0; i < n_mu(stock,year,season,age); i++) move(mu_row(cum_n_mu + i)-1,mu_col(cum_n_mu + i)-1) = one/(one + exp(-mu(mu_pointer(cum_n_mu + i)-1)));
      
      for(int i = 0; i < n_regions; i++) {
        vector<Type> trans_par = trans_mu.row(i);
        vector<int> do_move = can_move.row(i);
        vector<Type> pmove = additive_transform(trans_par, i, do_move);
        for(int j = 0; j < n_regions; j++) move(i,j) = pmove(j);
      }
  
}
  
