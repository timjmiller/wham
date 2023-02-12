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

template <class Type>
vector<matrix<Type> > get_selectivity(int n_years, int n_ages, int n_lengths, vector<Type> lengths, int n_selblocks, vector<matrix<Type> > selpars, vector<int> selblock_models)
{
  vector<matrix<Type> > selAL(n_selblocks);
  for(int b = 0; b < n_selblocks; b++)
  {
    matrix<Type> tmp(n_years, n_ages);
	matrix<Type> tmpL(n_years, n_lengths); // only for len selex
	vector<Type> newLengths(n_lengths); // to save length bin + 0.5*binsize
	Type binwidth = lengths(1) - lengths(0);
	tmp.setZero();
	tmpL.setZero();
    for(int l = 0; l < n_lengths; l++) newLengths(l) = lengths(l) + 0.5*binwidth; // saving length+0.5*binsize
    if(selblock_models(b) == 1) tmp = selpars(b); //proportions at age
    else
    { 
      if(selblock_models(b) == 2)
      { //increasing age-logistic
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
      { 
        if(selblock_models(b) == 3)
        {
		  // double logistic
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
          if(selblock_models(b) == 4)
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
		  } else { 
				if(selblock_models(b) == 5) { // age double normal
					  for(int y = 0; y < n_years; y++)
					  {
						Type p_1 = selpars(b)(y,0); // 
						Type p_2 = selpars(b)(y,1); // 
						Type p_3 = selpars(b)(y,2);
						Type p_4 = selpars(b)(y,3);
						Type p_5 = 1/(1+exp(-selpars(b)(y,4)));
						Type p_6 = 1/(1+exp(-selpars(b)(y,5)));
						Type age = 0.0;
						Type binwidth = 1;
						Type gammax = p_1 + binwidth + (0.99*n_ages - p_1 - binwidth)/(1 + exp(-p_2));
						Type alpha = 0.0;
						Type beta = 0.0;
						Type j_1 = 0.0;
						Type j_2 = 0.0;
						Type amin = 1;
						for (int a = 0; a < n_ages; a++)
						{
						  age += 1.0;
						  alpha = p_5 + (1 - p_5)*(exp(-pow(age - p_1, 2)/exp(p_3)) - exp(-pow(amin - p_1,2)/exp(p_3)))/(1-exp(-pow(amin - p_1,2)/exp(p_3)));
						  beta = 1 + (p_6 - 1)*(exp(-pow(age - gammax,2)/exp(p_4)) - 1)/(exp(-pow(n_ages - gammax,2)/exp(p_4)) - 1);
						  j_1 = 1/(1 + exp(-20*(age - p_1)/(1  + fabs(age - p_1))));
						  j_2 = 1/(1 + exp(-20*(age - gammax)/(1  + fabs(age - gammax))));
						  tmp(y,a) = alpha * (1 - j_1) + j_1*((1 - j_2) + j_2*beta);
						}
					  }
				} else {
					if(selblock_models(b) == 6) { // len-increasing logistic
						for(int y = 0; y < n_years; y++)
						{
						  Type l50 = selpars(b)(y,0); // l50 parameter in year y
						  Type k = selpars(b)(y,1); //  1/slope in year y
						  for(int l = 0; l < n_lengths; l++)
						  {
							tmpL(y,l) = 1.0/(1.0 + exp(-(newLengths(l) - l50)/k)); 
						  }
						  for(int l = 0; l < n_lengths; l++) tmpL(y,l) = tmpL(y,l)/tmpL(y,n_lengths-1); // standardize from 0 to 1?
						}
					} else { // len-decreasing logistic
						if(selblock_models(b) == 7) {
							for(int y = 0; y < n_years; y++)
							{
							  Type l50 = selpars(b)(y,0); // l50 parameter in year y
							  Type k = selpars(b)(y,1); //  1/slope in year y
							  for(int l = 0; l < n_lengths; l++)
							  {
								tmpL(y,l) = 1.0/(1.0 + exp((newLengths(l) - l50)/k)); 
							  }
							  for(int l = 0; l < n_lengths; l++) tmpL(y,l) = tmpL(y,l)/tmpL(y,n_lengths-1); // standardize from 0 to 1?
							}
						} else { // length double normal
						  for(int y = 0; y < n_years; y++)
						  {				
							Type p_1 = selpars(b)(y,0); // 
							Type p_2 = selpars(b)(y,1); // 
							Type p_3 = selpars(b)(y,2);
							Type p_4 = selpars(b)(y,3);
							Type p_5 = 1/(1+exp(-selpars(b)(y,4)));
							Type p_6 = 1/(1+exp(-selpars(b)(y,5)));
							Type lmax = max(newLengths);
							Type lmin = min(newLengths);
							Type gammax = p_1 + binwidth + (0.99*lmax - p_1 - binwidth)/(1 + exp(-p_2));
							Type alpha = 0.0;
							Type beta = 0.0;
							Type j_1 = 0.0;
							Type j_2 = 0.0;
							for (int l = 0; l < n_lengths; l++)
							{
							  alpha = p_5 + (1 - p_5)*(exp(-pow(newLengths(l) - p_1, 2)/exp(p_3)) - exp(-pow(lmin - p_1,2)/exp(p_3)))/(1-exp(-pow(lmin - p_1,2)/exp(p_3)));
							  beta = 1 + (p_6 - 1)*(exp(-pow(newLengths(l) - gammax,2)/exp(p_4)) - 1)/(exp(-pow(lmax - gammax,2)/exp(p_4)) - 1);
							  j_1 = 1/(1 + exp(-20*(newLengths(l) - p_1)/(1  + fabs(newLengths(l) - p_1))));
							  j_2 = 1/(1 + exp(-20*(newLengths(l) - gammax)/(1  + fabs(newLengths(l) - gammax))));
							  tmpL(y,l) = alpha * (1 - j_1) + j_1*((1 - j_2) + j_2*beta);
							}
						  }	
						}
					}
					
				}
		  }
        }
      }
    }
	
    if(selblock_models(b) < 6) selAL(b) = tmp; // age selex
	else selAL(b) = tmpL; // len selex
	
  }
  return selAL;
}

template <class Type>
vector<Type> get_waa_y(array<Type> waa, int y, int na, int pointer){
  vector<Type> waay(na);
  for(int a = 0; a < na; a++) waay(a) = waa(pointer-1, y, a);
  return(waay);
}

template <class Type>
matrix<Type> get_waacatch_y(array<Type> waa, int y, int na, vector<int> pointer){
  matrix<Type> waay(pointer.size(), na);
  for(int f = 0; f< pointer.size(); f++) for(int a = 0; a < na; a++) waay(f,a) = waa(pointer(f)-1, y, a);
  return(waay);
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
Type get_SPR(Type log_F, vector<Type> M, matrix<Type> sel, vector<Type> mat, vector<Type> waassb, Type fracyearSSB)
{
  int n_ages = M.size();
  int n_fleets = sel.rows();
  Type SPR = 0.0, ntemp = 1.0;
  vector<Type> Z = M;
  matrix<Type> F = exp(log_F) * sel;
  //Z = F + M;
  for(int age=0; age<n_ages-1; age++)
  {
    for(int f = 0; f < n_fleets;f++) Z(age) += F(f,age);
    SPR += ntemp * mat(age) * waassb(age) * exp(-fracyearSSB * Z(age));
    ntemp *= exp(-Z(age));
  }
  for(int f = 0; f < n_fleets;f++) Z(n_ages-1) += F(f,n_ages-1);
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
Type get_YPR(Type log_F, vector<Type> M, matrix<Type> sel, matrix<Type> waacatch)
{
  int n_ages = M.size();
  int n_fleets = sel.rows();
  Type YPR = 0.0, ntemp = 1.0;
  vector<Type> Z = M;

  matrix<Type> F = exp(log_F) * sel;
  for(int age=0; age<n_ages; age++) for(int f = 0; f < n_fleets; f ++) Z(age) += F(f,age);
  for(int age=0; age<n_ages-1; age++) {
    for(int f = 0; f < n_fleets; f ++) YPR += ntemp * F(f,age) * waacatch(f,age) * (1.0 - exp(-Z(age)))/Z(age);
    ntemp *= exp(-Z(age));
  }
  ntemp /= 1.0 - exp(-Z(n_ages-1));
  for(int f = 0; f < n_fleets; f ++) YPR += ntemp * F(f,n_ages-1) * waacatch(f,n_ages-1) * (1.0 - exp(-Z(n_ages-1)))/Z(n_ages-1);
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

/* calculate beverton-holt or ricker equilibrium yield */
/* uses selectivity and waa by fleet */
template<class Type>
struct sr_yield_fleet {
  /* Data and parameter objects for yield calculation: */
  Type SR_a;
  Type SR_b;
  vector<Type> M;
  matrix<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  matrix<Type> waacatch;
  Type fracyearSSB;
  int sr_type;

  /* Constructor */
  sr_yield_fleet(Type SR_a_, Type SR_b_,
  vector<Type> M_,
  matrix<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  matrix<Type> waacatch_,
  Type fracyearSSB_, int sr_type_) :
    SR_a(SR_a_), SR_b(SR_b_), M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), waacatch(waacatch_), fracyearSSB(fracyearSSB_), sr_type(sr_type_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    int n_fleets = sel.rows();
    T YPR = 0, SPR = 0, ntemp = 1, R;
    vector<T> Z = M.template cast<T>();

    matrix<T> F = exp(log_F(0)) * sel.template cast<T>();
    //Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++) for(int f = 0; f < n_fleets; f++) Z(age) += F(f,age);
    for(int age=0; age<n_ages-1; age++) {
      for(int f = 0; f < n_fleets; f++) YPR += ntemp * F(f,age) * T(waacatch(f,age)) * (1- exp(-Z(age)))/Z(age);
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    for(int f = 0; f < n_fleets; f++) YPR += ntemp * F(f,n_ages-1) * T(waacatch(f,n_ages-1)) * (1 - exp(-Z(n_ages-1)))/Z(n_ages-1);
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

/* calculate SPR at F */
/* uses selectivity and by fleet */
template<class Type>
struct spr_F_fleet {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  matrix<Type> sel;
  vector<Type> mat;
  vector<Type> waassb;
  Type fracyearSSB;

  /* Constructor */
  spr_F_fleet(
  vector<Type> M_,
  matrix<Type> sel_,
  vector<Type> mat_,
  vector<Type> waassb_,
  Type fracyearSSB_) :
    M(M_), sel(sel_), mat(mat_),
  waassb(waassb_), fracyearSSB(fracyearSSB_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it maximizes yield
    int n_ages = M.size();
    int n_fleets = sel.rows();
    T SPR = 0, ntemp = 1;
    vector<T> Z = M.template cast<T>();

    matrix<T> F = exp(log_F(0)) * sel.template cast<T>();
    for(int age=0; age<n_ages; age++) for(int f = 0; f < n_fleets; f++) Z(age) += F(f,age);
    //Z = F + M.template cast<T>();
    for(int age=0; age<n_ages-1; age++) {
      SPR += ntemp * T(mat(age) * waassb(age)) * exp(-T(fracyearSSB) * Z(age));
      ntemp *= exp(-Z(age));
    }
    ntemp /= 1 - exp(-Z(n_ages-1));
    SPR += ntemp * T(mat(n_ages-1) * waassb(n_ages-1)) * exp(-T(fracyearSSB)*Z(n_ages-1));

    return SPR;
  }
};

template <class Type>
matrix<Type> get_SPR_res(matrix<Type> MAA, array<Type> FAA, vector<int> which_F_age, array<Type> waa, int waa_pointer_ssb, vector<int> waa_pointer_fleets,
  matrix<Type> mature, Type percentSPR, vector<Type> predR, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<Type> F_init)
{
  int n = 10;
  int ny = MAA.rows();
  int na = MAA.cols();
  int nf = waa_pointer_fleets.size();
  matrix<Type> res(ny, 5+n);
  vector<Type> log_FXSPR(ny), log_FXSPR_i(1), log_SPR_FXSPR(ny);
  vector<Type> log_SSB_FXSPR(ny), log_Y_FXSPR(ny);
  matrix<Type> log_FXSPR_iter(ny,n);
  //matrix<Type> catchWAA = extract_matrix_array3(waa, waa_pointer_tot_catch-1);
  matrix<Type> ssbWAA = extract_matrix_array3(waa, waa_pointer_ssb-1);
  vector<Type> log_YPR_FXSPR(ny), waassb(na), mat(na), M(na);
  matrix<Type> sel(nf,na), waacatch(nf, na);
  for(int y = 0; y < ny; y++)
  {
    log_FXSPR_iter(y,0) = res(y,5) = log(F_init(y));
    M = MAA.row(y);
    waassb = ssbWAA.row(y);
    waacatch = get_waacatch_y(waa, y, na, waa_pointer_fleets);
    //waacatch = catchWAA.row(y);
    mat = mature.row(y);
    Type F_full = 0;
    for(int f = 0; f < nf; f++) F_full += FAA(y,f,which_F_age(y)-1);
    for(int f = 0; f < nf; f++) for(int a = 0; a < na; a++) sel(f,a) = FAA(y,f,a)/F_full;
    spr_F_fleet<Type> sprF(M, sel, mat, waassb, fracyr_SSB(y));
    for (int i=0; i<n-1; i++)
    {
      log_FXSPR_i(0) = log_FXSPR_iter(y,i);
      vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
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
    res(y,3) = log_SPR_FXSPR(y);
    res(y,4) = log_YPR_FXSPR(y);
  }
  return res;
}

template <class Type>
vector<Type> get_static_SPR_res(matrix<Type> MAA, array<Type> FAA, int which_F_age, array<Type> waa, int waa_pointer_ssb, vector<int> waa_pointer_fleets,
  matrix<Type> mature, Type percentSPR, matrix<Type> NAA, vector<Type> fracyr_SSB, Type F_init, 
  vector<int> years_M, vector<int> years_mat, vector<int> years_sel, vector<int> years_waa_ssb, vector<int> years_waa_catch, vector<int> years_R)
{
  int n = 10;
  int na = MAA.cols();
  int nf = waa_pointer_fleets.size();
  //matrix<Type> waacatch = extract_matrix_array3(waa, waa_pointer_tot_catch-1);
  matrix<Type> waassb = extract_matrix_array3(waa, waa_pointer_ssb-1);
  matrix<Type> sel(nf,na), waa_c(nf,na);
  vector<Type> M(na), mat(na), waa_s(na);
  sel.setZero(); M.setZero(); mat.setZero(); waa_s.setZero(), waa_c.setZero();
  Type ssbfrac = 0, R = 0;
  //get average inputs over specified years
  for(int y = 0; y < years_R.size(); y++) R += NAA(years_R(y),0)/years_R.size();
  for(int y = 0; y < years_mat.size(); y++) for(int a = 0; a < na; a++) mat(a) += mature(years_mat(y),a)/years_mat.size();
  for(int y = 0; y < years_M.size(); y++) for(int a = 0; a < na; a++) M(a) += MAA(years_M(y),a)/years_M.size();
  for(int y = 0; y < years_waa_catch.size(); y++) for(int f = 0; f < nf; f++) for(int a = 0; a < na; a++) {
    waa_c(f,a) += waa(waa_pointer_fleets(f)-1, years_waa_catch(y),a)/years_waa_catch.size();
  }
  for(int y = 0; y < years_waa_ssb.size(); y++) {
   ssbfrac += fracyr_SSB(years_waa_ssb(y))/years_waa_ssb.size();
   for(int a = 0; a < na; a++) waa_s(a) += waassb(years_waa_ssb(y),a)/years_waa_ssb.size();
  }
  for(int y = 0; y < years_sel.size(); y++) for(int a = 0; a < na; a++) for(int f = 0; f < nf; f++) {
    sel(f,a) += FAA(years_sel(y),f,a)/years_sel.size();
  }
  Type F_mult= 0;
  for(int f = 0; f < nf; f++) F_mult += sel(f,which_F_age-1);
  sel = sel/F_mult;
  Type spr0 = get_SPR_0(M, mat, waa_s, ssbfrac); 

  vector<Type> res(6+n), log_FXSPR_i(1), log_FXSPR_iter(n);
  log_FXSPR_iter(0) = res(6) = log(F_init);
  spr_F_fleet<Type> sprF(M, sel, mat, waa_s, ssbfrac);
  for (int i=0; i<n-1; i++)
  {
    log_FXSPR_i(0) = log_FXSPR_iter(i);
    vector<Type> grad_spr_F = autodiff::gradient(sprF,log_FXSPR_i);
    log_FXSPR_iter(i+1) = log_FXSPR_iter(i) - (sprF(log_FXSPR_i) - 0.01*percentSPR * spr0)/grad_spr_F(0);// /hess_sr_yield(0,0);
    res(6+i+1) = log_FXSPR_iter(i+1);
  }
  res(0) = log_FXSPR_iter(n-1); //log_FXSPR
  res(3) = log(get_SPR(res(0), M, sel, mat, waa_s, ssbfrac)); //log_SPR_FXSPR, should be log(X*SPR0/100)
  res(1) = log(R) + res(3); //log_SSB_FXSPR
  res(4) = log(get_YPR(res(0), M, sel, waa_c)); //log_YPR_FXSPR
  res(2) = log(R) + res(4); //log_Y_FXSPR
  res(5) = log(spr0); //log_SPR0
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

/* calculate total catch at F across multiple fleets*/
template<class Type>
struct log_catch_fleets_F {
  /* Data and parameter objects for calculation: */
  vector<Type> M;
  matrix<Type> sel; //n_fleets x n_ages
  matrix<Type> waacatch; //n_fleets x n_ages
  vector<Type> NAA;

  /* Constructor */
  log_catch_fleets_F(
  vector<Type> M_,
  matrix<Type> sel_,
  matrix<Type> waacatch_,
  vector<Type> NAA_) :
    M(M_), sel(sel_), waacatch(waacatch_), NAA(NAA_) {}

  template <typename T> //I think this allows you to differentiate the function wrt whatever is after operator() on line below
  T operator()(vector<T> log_F) { //find such that it achieves required catch
    int n_ages = M.size();
    int n_fleets = waacatch.rows();
    T Catch = 0;
    vector<T> Z = M.template cast<T>();

    matrix<T> F = exp(log_F(0)) * sel.template cast<T>();
    //Z = F + M.template cast<T>();
    for(int age=0; age<n_ages; age++) {
      for(int f = 0; f < n_fleets; f++) Z(age) += F(f,age);
      for(int f = 0; f < n_fleets; f++){
        Catch += T(waacatch(f,age)) * T(NAA(age)) * F(f,age) *(1 - exp(-Z(age)))/Z(age);
      }
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

//multiple fleets
template <class Type>
Type get_F_from_log_Catch(Type Catch, vector<Type> NAA, vector<Type> M, matrix<Type> sel, matrix<Type> waacatch, Type F_init)
{
  int n = 10;
  vector<Type> log_F_i(1);
  vector<Type> log_F_iter(n);
  log_F_iter.fill(log(F_init)); //starting value
  log_catch_fleets_F<Type> logcatchF(M, sel, waacatch, NAA);
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

//fleet-specific FAA
template <class Type>
Type get_FMSY(Type log_a, Type log_b, vector<Type> M, matrix<Type> sel, matrix<Type> waacatch, vector<Type> waassb,
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
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: matrix of 0/1 with first column determining if Ecov is affecting recruitment
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: array that holds linear predictor for Ecov
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
      int nE = Ecov_where.rows();
      for(int i=0; i < nE; i++){
        if(Ecov_where(i,0) == 1) if(Ecov_how(i) == 1) pred_recruit *= exp(Ecov_lm(i,0,y,0));
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
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm, matrix<Type> ZAA){

  /*
   * y: year (between 1 and n_years_model+n_years_proj)
   * recruit_model: which recruitment model (1-4)
   * mean_rec_pars: vector of any recruitment parameters (defined in main code)
   * SSB: vector of yearly SSB (uses y-1 for any S-R relationship)
   * NAA: matrix of numbers at age
   * log_SR_a: yearly "a" parameters for SR function
   * log_SR_b: yearly "b" parameters for SR function
   * Ecov_where: matrix of 0/1 integer with first column determining if Ecov i affects recruitment (= 1)
   * Ecov_how: integer vector with an element that tells how the Ecov is affecting recruitment
   * Ecov_lm: array that holds linear predictor for Ecov
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
Type get_SSB(matrix<Type> NAA, matrix<Type> ZAA, array<Type> waa, matrix<Type> mature, 
			 //matrix<Type> mature_len, int weight_model, array<Type> out_phi_mat, matrix<Type> watl, 
			 int y, int ssbpointer, vector<Type> fracyr_SSB){
	Type SSB = 0;
	//Type fec_age = 0;
	//int n_lengths = mature_len.cols();
	int na = mature.cols();
    // Calculate SSB using mature (at age) or mature_len
	//if((weight_model == 1) | (weight_model == 3)) { // use age information
	  vector<Type> waay = get_waa_y(waa,y,na,ssbpointer);
	  for(int a = 0; a < na; a++) {
		  SSB += NAA(y,a) * waay(a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
	  }
	//}
	//if(weight_model == 2) { // use len information
	//	for(int a = 0; a < na; a++) {
	//		fec_age = 0;
	//		for(int l = 0; l < n_lengths; l++) fec_age += out_phi_mat(y,l,a)*watl(y,l)*mature_len(y,l)*mature(y,a);
	//		SSB += NAA(y,a) * fec_age * exp(-ZAA(y,a)*fracyr_SSB(y));
	//	}
	//}
  return(SSB);
}
 
template <class Type>
matrix<Type> get_F_proj(int y, int n_fleets, vector<int> proj_F_opt, array<Type> FAA, matrix<Type> NAA, matrix<Type> MAA, matrix<Type> mature, 
  matrix<Type> waacatch, vector<Type> waassb, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, vector<int> which_F_age, Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, Type F_init, 
  vector<Type> log_a, vector<Type> log_b, int recruit_model, Type percentFMSY){
    /* 
     get F to project for next time step
              y:  year of projection (>n_years_model)
     proj_F_opt:  for each year, how to specify F for projection (1 to 6)
            FAA:  FAA array from main code.
            NAA:  NAA matrix from main code
            MAA:  MAA matrix from main code.
         mature:  maturity matrix from main code.
       waacatch:  matrix (n_fleets x n_ages) weight at age in catch to use for year y (use function get_waacatch_y defined above)
         waassb:  weight at age for SSB to use in year y (use function get_waa_y defined above)
     fracyr_SSB:  vector of yearly fractions of the year when spawning occurs
       log_SPR0:  vector of yearly log(unfished SSB/R) 
  avg_years_ind:  integer vector of years to average F for projection
  n_years_model:  number of years before projection begins
    which_F_age:  define which age has max F
     percentSPR:  percentage (0-100) of unfished spawning potential to determine F_percentSPR
    proj_Fcatch:  vector (n_years_proj) of user specified Fishing mortality rates to project
   percentFXSPR:  percentage (0-100) of F_percentSPR to use in catch, e.g. GOM cod uses F = 75% F_40%SPR
         F_init:  initial value to use for FXSPR or FMSY newton method
          log_a:  annual log(a) for stock-recruit relationship
          log_b:  annual log(b) for stock-recruit relationship
  recruit_model:  integer for which type of recruit model is assumed (= 3 or 4 for using Fmsy)
    percentFMSY:  percentage (0-100) of FMSY to use in catch.
    */
  int n_toavg = avg_years_ind.size();
  int n_ages = waacatch.cols();
  int proj_F_opt_y = proj_F_opt(y-n_years_model);

  //proj_F_opt == 1, last year F (default)
  matrix<Type> FAA_proj(n_fleets, n_ages);
  vector<Type> FAA_tot_proj(n_ages);
  FAA_tot_proj.setZero();
  FAA_proj.setZero();
  if(proj_F_opt_y == 1){ // last year F (default)
    for(int f = 0; f < n_fleets; f++) for(int a = 0; a < n_ages; a++) FAA_proj(f,a) = FAA(n_years_model-1,f,a);
  }
  else { //proj_F_opt_y>1
    //option 1: average F is by fleet and Ftot is sum of fleet averages
    //when there is more than 1 fleet, the sum of predicted catch across fleets will not generally equal the total catch using FAA_tot and waa_totcatch.
    for(int f = 0; f < n_fleets; f++) 
    {
      for(int a = 0; a < n_ages; a++) for(int i = 0; i < n_toavg; i++){
        FAA_proj(f,a) += FAA(avg_years_ind(i),f,a);
      }
      FAA_proj.row(f) /= Type(n_toavg);
    }
      
    //proj_F_opt_y == 2, average F
    if(proj_F_opt_y == 2){
      //already defined
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
    Type fracyr_SSB_y = fracyr_SSB(y);
    Type log_SPR0_y = log_SPR0(y);
    Type log_a_y = log_a(y);
    Type log_b_y = log_b(y);
    if(proj_F_opt_y == 3){ // F at X% SPR
      FAA_tot_proj = get_FXSPR(M, sel_tot_proj, waassb, mat, percentSPR, fracyr_SSB_y, log_SPR0_y, F_init) * sel_tot_proj * 0.01*percentFXSPR;
      FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
    }
    if(proj_F_opt_y == 4){ // user-specified F
      if(proj_Fcatch(y-n_years_model) < 1e-10){ // if F = 0, sel_proj is NaN
        FAA_proj.setZero();
      } else {
        FAA_proj = Type(proj_Fcatch(y-n_years_model)) * sel_proj;
        FAA_tot_proj = FAA_proj.colwise().sum();
        //FAA_tot_proj = Type(proj_Fcatch(y-n_years_model)) * sel_tot_proj;
        //FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
      }
    }
    if(proj_F_opt_y == 5){ // calculate F from user-specified catch
      Type thecatch = proj_Fcatch(y-n_years_model);
      if(thecatch < 1e-10){ // if catch = 0, F = 0 and sel_proj is NaN
        FAA_proj.setZero();
      } else {
        //FAA_tot_proj = get_F_from_log_Catch(thecatch, NAA_y, M, sel_tot_proj, waacatch, F_init) * sel_tot_proj;
        FAA_proj = get_F_from_log_Catch(thecatch, NAA_y, M, sel_proj, waacatch, F_init) * sel_proj;
        FAA_tot_proj = FAA_proj.colwise().sum();
      }
    }
    if(proj_F_opt_y == 6){ //Fmsy
      FAA_proj = get_FMSY(log_a_y, log_b_y, M, sel_proj, waacatch, waassb, mat, fracyr_SSB_y, log_SPR0_y, recruit_model, F_init) * sel_proj * 0.01* percentFMSY;
      FAA_tot_proj = FAA_proj.colwise().sum();
      //FAA_proj = FAA_tot_proj(which_F_age(y)-1) * sel_proj;
    }
                             
  }
  return(FAA_proj);
}


template <class Type>
matrix<Type> sim_pop(array<Type> NAA_devs, int recruit_model, vector<Type> mean_rec_pars, vector<Type> SSBin, matrix<Type> NAAin, vector<Type> log_SR_a, 
  vector<Type> log_SR_b, matrix<int> Ecov_where, vector<int> Ecov_how, array<Type> Ecov_lm, int n_NAA_sigma, 
  int do_proj, vector<int> proj_F_opt, array<Type> FAA, matrix<Type> FAA_tot, matrix<Type> MAA, matrix<Type> mature, array<Type> waa, 
  vector<int> waa_pointer_fleets, int waa_pointer_ssb, vector<Type> fracyr_SSB, vector<Type> log_SPR0, vector<int> avg_years_ind, 
  int n_years_model, int n_fleets, vector<int> which_F_age, Type percentSPR, vector<Type> proj_Fcatch, Type percentFXSPR, vector<Type> F_proj_init, Type percentFMSY){

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
        //vector<Type> waacatch = get_waa_y(waa, y, n_ages, waa_pointer_totcatch);
        matrix<Type> waacatch = get_waacatch_y(waa, y, n_ages, waa_pointer_fleets);
        vector<Type> waassb = get_waa_y(waa, y, n_ages, waa_pointer_ssb);
        matrix<Type> FAA_proj = get_F_proj(y, n_fleets, proj_F_opt, FAA, NAA, MAA, mature, waacatch, waassb, fracyr_SSB, log_SPR0, avg_years_ind, n_years_model,
         which_F_age, percentSPR, proj_Fcatch, percentFXSPR, F_proj_init(y-n_years_model), log_SR_a, log_SR_b, recruit_model, percentFMSY);
        for(int f = 0; f < n_fleets; f++) for(int a = 0; a< n_ages; a++) FAA(y,f,a) = FAA_proj(f,a);
        FAA_tot.row(y) = FAA_proj.colwise().sum();
        ZAA.row(y) = FAA_tot.row(y) + MAA.row(y);
      }
    } // end proj F
    SSB(y) = get_SSB(NAA,ZAA,waa,mature,y,waa_pointer_ssb,fracyr_SSB);
  } // end pop model loop
  
  matrix<Type> out(ny, 2 * n_ages + n_fleets * n_ages + 1);
  for(int i = 0; i < n_ages; i++) out.col(i) = NAA.col(i);
  for(int i = n_ages; i < 2 * n_ages; i++) out.col(i) = pred_NAA.col(i-n_ages);
  for(int y = 0; y < ny; y++) for(int i = 2 ; i < (2 + n_fleets); i++) for(int a = 0; a < n_ages; a++) out(y,i*n_ages + a) = FAA(y,i-2,a);
  out.col(out.cols()-1) = SSB;
  return(out);
}

template <class Type>
array<Type> pred_LAA(matrix<Type> mLAA_jan1, int n_yrs, int n_years_model, vector<Type> expSD,
					  array<Type> GW_par, vector<Type> lengths, vector<Type> fracyr_vec, int growth_model, Type age_L1){

  Type Lminp = min(lengths);
  Type Lmaxp = max(lengths);
  Type Fac1 = 0.0;
  Type Fac2 = 0.0;
  Type Ll1p = 0.0;
  Type Llp = 0.0;
  int n_ages = mLAA_jan1.cols();
  int n_lengths = lengths.size();
  array<Type> out(n_yrs, n_lengths, n_ages);
  vector<Type> mLAA(n_ages);
  Type Grate = 0.0;
  Type b_len = 0.0;
  Type last_linear = 0.0;
  Type this_SD = 0.0;
  Type Slope = 0.0;

	for(int y = 0; y < n_yrs; y++) {
		  mLAA.setZero();
		  int y_1 = y + 1;
		  int yuse = y; // used in fraction information
		  if(y > (n_years_model - 1)) yuse = n_years_model -1; //some things only go up to n_years_model-1
		  if(y == (n_yrs - 1)) y_1 = y;
		  Type fracyr = fracyr_vec(yuse);
		  for(int a = 0; a < n_ages; a++)
		  {  
			if(growth_model == 1) { // classic von Bertalanffy curve
				b_len = (GW_par(y,a,2) - Lminp)/age_L1;
				if((a + 1.0 + fracyr) <= age_L1) { // use linear growth
					mLAA(a) = Lminp + b_len*(a+1.0+fracyr);  
				} else { // other options
					if((a+1.0) >= age_L1) { // only growth curve
						if(a == (n_ages-1)) mLAA(a) = mLAA_jan1(y,a); // no growth for age plus group
						else mLAA(a) = mLAA_jan1(y,a) + (mLAA_jan1(y,a) - GW_par(y,a,1))*(exp(-GW_par(y,a,0)*fracyr) - 1.0);
					} else { // linear + growth curve mixed
						last_linear = Lminp + b_len*age_L1; // last age (cont) with linear growth 
						mLAA(a) = last_linear + (last_linear - GW_par(y,a,1))*(exp(-GW_par(y,a,0)*(a+1.0+fracyr-age_L1)) - 1.0); 
					}
				}
			}
			if(growth_model == 2) { // Richards growth curve	
				b_len = (GW_par(y,a,2) - Lminp)/age_L1;			
				if((a + 1.0 + fracyr) <= age_L1) { // use linear growth
					mLAA(a) = Lminp + b_len*(a+1.0+fracyr);  
				} else { // other options
					if((a+1.0) >= age_L1) { // only growth curve
						if(a == (n_ages-1)) mLAA(a) = mLAA_jan1(y,a); // no growth for age plus group
						else mLAA(a) = pow(pow(mLAA_jan1(y,a),GW_par(y,a,3)) + (pow(mLAA_jan1(y,a),GW_par(y,a,3)) - pow(GW_par(y,a,1),GW_par(y,a,3)))*(exp(-GW_par(y,a,0)*fracyr) - 1.0),1/GW_par(y,a,3));
					} else { // linear + growth curve mixed
						last_linear = Lminp + b_len*age_L1; // last age (cont) with linear growth 
						mLAA(a) = pow(pow(last_linear,GW_par(y,a,3)) + (pow(last_linear,GW_par(y,a,3)) - pow(GW_par(y,a,1),GW_par(y,a,3)))*(exp(-GW_par(y,a,0)*(a+1.0+fracyr-age_L1)) - 1.0),1/GW_par(y,a,3));
					}
				}
			}
			if(growth_model == 3){ // nonparametric approach
				if(a < (n_ages-1)) { // for a < n_ages - 1, linear interpolation
					Grate = (mLAA_jan1(y_1,a+1) - mLAA_jan1(y,a))*fracyr;
					mLAA(a) = mLAA_jan1(y,a) + Grate;	
				} else { // last age. use last Grate
					mLAA(a) = mLAA_jan1(y,a) + Grate;
				}
			}
			
			// SD calculation (again): 
			if(growth_model < 3) { // for parametric approach
				if((a + 1.0 + fracyr) <= age_L1) { // same as SD1
					this_SD = expSD(0); 
				} else { 
					if(a == (n_ages-1)) { // same as SDA
						this_SD = expSD(1);
					} else { // linear interpolation
						Slope = (expSD(1) - expSD(0))/(GW_par(y,a,1)-GW_par(y,a,2));
						this_SD = expSD(0) + Slope*(mLAA(a)-GW_par(y,a,2));  
					}
				}				
			}
			if(growth_model == 3) { // only works for year effect
				Slope = (expSD(1) - expSD(0))/(mLAA_jan1(y,n_ages-1)-mLAA_jan1(y,0));
				this_SD = expSD(0) + Slope*(mLAA(a)-mLAA_jan1(y,0));  
			}

			for(int l = 0; l < n_lengths; l++) {
				
				if(l == 0) { 
					Fac1 = (Lminp - mLAA(a))/this_SD;
					out(y,l,a) = pnorm(Fac1);  
				} else {
					if(l == (n_lengths-1)) { 
						Fac1 = (Lmaxp - mLAA(a))/this_SD;
						out(y,l,a) = 1.0 - pnorm(Fac1);  
					} else { 
						Ll1p = lengths(l+1);
						Llp = lengths(l);
						Fac1 = (Ll1p - mLAA(a))/this_SD;
						Fac2 = (Llp - mLAA(a))/this_SD;
						out(y,l,a) = pnorm(Fac1) - pnorm(Fac2);  
					}
				}
			}
		  }
	}
	
  return(out);
}

template <class Type>
matrix<Type> get_fracyr_WAA(vector<Type> WAA_jan1, vector<Type> WAA_jan1_y1, Type fracyr){
  int n_ages = WAA_jan1.size();
  vector<Type> WAA(n_ages);
  Type Grate = 0.0;
  
  	for(int a = 0; a < n_ages; a++)
	 {  
		if(a < (n_ages-1)) { // for a < n_ages - 1. See Crane et al 2019
			Grate = pow(WAA_jan1_y1(a+1)/WAA_jan1(a), 1.0/365.0) - 1.0;
			WAA(a) = WAA_jan1(a)*pow(1.0 + Grate, fracyr*365.0);
		} else { // for last age, use last Grate
			WAA(a) = WAA_jan1(a)*pow(1.0 + Grate, fracyr*365.0); 
		}
		
	}
	  
	return(WAA);
}

// Get selectivity at age from selectivity at length
// Only used when selectivity at length is specified
// Uses phi_matrix at fracyr
template <class Type>
vector<Type> get_selAA_from_selAL(matrix<Type> selAL, int y, int this_sel_model, array<Type> fracyr_phi_mat){
	  
	  vector<int> dims = fracyr_phi_mat.dim;
	  int n_ages = dims(2);
	  int n_lengths = dims(1);
	  vector<Type> selAA(n_ages); // n_ages
	  if(this_sel_model < 6) { // for age models
		  selAA = selAL.row(y);  // same as calculated in selAL
	  } else { // transform for all years
		  for(int a = 0; a < n_ages; a++) {
		  Type sumSelex = 0.0;
			  for(int l = 0; l < n_lengths; l++) {
				sumSelex += fracyr_phi_mat(y,l,a)*selAL(y,l);
			  }
		  selAA(a) = sumSelex;
		  }
	  }
	return(selAA);
}
