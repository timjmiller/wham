

/* calculate SPR at F */
template<class Type>
struct spr_F_sp {
  /* Data and parameter objects for calculation: */
  array<Type> M;
  array<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  Type fracyearSSB;

  /* Constructor */
  spr_F_sp(
  matrix<Type> M_,
  matrix<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  Type fracyearSSB_) :
    M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), fracyearSSB(fracyearSSB_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    T SPR = 0, ntemp = 1;
    array<T> F = exp(log_F(0)) * sel.template cast<T>();
    //array<T> Z = F + M.template cast<T>();
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
