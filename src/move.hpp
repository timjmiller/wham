nll_mu_prior_re

simulate_mu_prior_re
  
nll_mu_re
  

simulate_mu_re
  
array<Type> get_mu(trans_mu, mu_re, can_move, mig_type)
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
  
