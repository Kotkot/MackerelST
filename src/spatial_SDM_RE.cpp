////////////////////////////////////////////////////
// Description:
// Model for the spatial CPUE by age i.e CPUE(age) at a given location
// Each age class follow a Tweedie distribution (or we can decide per age group) 
// AR1 process for between age correlation (distance is calculated as age difference)
// Each age group has its own spatial field
// we can model the covariate effect separaetely for each group

// Author: 
// Kotaro Ono
//
// Version: 1
// Detail: 
//
// Important note: 
// 1. No barrier effect by year
//


#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Function to important barrier-SPDE code
template<class Type>
struct spde_barrier_t{
  vector<Type> C0;
  vector<Type> C1;
  Eigen::SparseMatrix<Type> D0;
  Eigen::SparseMatrix<Type> D1;
  Eigen::SparseMatrix<Type> I;
  spde_barrier_t(SEXP x){           // x = List passed from R 
    C0 = asVector<Type>(getListElement(x,"C0"));
    C1 = asVector<Type>(getListElement(x,"C1"));
    D0 = tmbutils::asSparseMatrix<Type>(getListElement(x,"D0"));
    D1 = tmbutils::asSparseMatrix<Type>(getListElement(x,"D1"));
    I = tmbutils::asSparseMatrix<Type>(getListElement(x,"I"));
  }
};

// Function to calculate Q (precision) matrix using barrier-SPDE
template<class Type>
Eigen::SparseMatrix<Type> Q_spde(spde_barrier_t<Type> spde_barrier, Type kappa, vector<Type> c){
  //using namespace Eigen;
  vector <Type> range(2);
  range(0) = sqrt(8.0)/kappa*c(0);
  range(1) = range(0)*c(1);
  Type pi = 3.141592;
  
  int dimLatent = spde_barrier.D0.row(0).size();
  vector<Type> Cdiag(dimLatent);
  Eigen::SparseMatrix<Type > Cinv(dimLatent,dimLatent);
  
  Cdiag = spde_barrier.C0*pow(range(0),2.0) + spde_barrier.C1*pow(range(1),2.0);
  for(int i =0; i<dimLatent; ++i){
    Cinv.coeffRef(i,i) = 1/Cdiag(i);
  }
  
  Eigen::SparseMatrix<Type>A = spde_barrier.I;
  A = A + (pow(range(0),2.0)/8.0) * spde_barrier.D0 + (pow(range(1),2.0)/8.0) * spde_barrier.D1;
  
  Eigen::SparseMatrix<Type> Q = A.transpose() * Cinv * A/pi *2 * 3;
  
  return Q;
}

// Repeat vector  
template <class Type>
vector<Type> RepeatVector(vector<Type> x, int times)
{
  int n = x.size() * times;
  vector<Type> res(n);
  int k = 0;
  for (int i = 0; i < times; i++) {
    for (int j = 0; j < x.size(); j++) {
      res[k] = x(j);
      k++;
    }
  }
  return res;
}

// Parameter transform for the autocorrelation coefficient
// approach 1
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(0.5) * x)) - Type(1);}

// approach 2
template <class Type>
Type zerofive_to_one(Type x)
{
  return Type(1) - Type(0.5) * invlogit(x);
}


// some likelihood functions
template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0)
{
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0)
{
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type logistic_func(Type x, Type s50, Type slope) {
  // logistic function controling the density-dependent effect. 
  // logistic function is similar to length or size based selectivity
  // in fisheries, parameterized by the points at which f(x) = 0.5 and the width of the slope
  Type pred = Type(1.0) / (Type(1.0) + exp(-log(Type(19.0)) * (x - s50) / slope));
  return pred;
}


// some specifications of available options
enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  inverse_link  = 3,
  boxcox_link   = 4
};

template <class Type>
Type InverseLink(Type eta, int link, Type powcoef)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = exp(eta);
    break;
  case logit_link:
    out = invlogit(eta);
    break;
  case inverse_link:
    out = Type(1.0) / eta;
    break;
  case boxcox_link:
    if (powcoef == Type(0)) out = exp(eta);
    if (powcoef != Type(0)) out = pow(powcoef * eta + 1, 1/powcoef);
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

enum valid_family {
  tweedie_family  = 0
  // gaussian_family   = 1,
  // poisson_family    = 2,
  // gamma_family      = 3,
  // nb_family         = 4,
  // student_family    = 5,
  // lognormal_family  = 6
};



template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  
  // Data section
  DATA_INTEGER(Nage);                         //The number of age groups
  DATA_MATRIX(X);                             //Design matrix for the fixed effects (same for each age groups)
  // DATA_VECTOR(yobs_vec);                   //the observed age samples at each station  
  // DATA_IMATRIX(yobs_index);                // indexing the vector of observation for the above one-step predictions;
  DATA_ARRAY(yobs);                           //the observed age samples at each station
  DATA_INTEGER(Nobs);                         //Number of observations
  DATA_INTEGER(spatialmodel);                 //whether it is a spatial model or not
  // Barrier effect
  DATA_STRUCT(spde_barrier,spde_barrier_t);   		//Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  DATA_VECTOR(barrier_scaling);               //scaling of range

  // Random intercepts:
  DATA_IMATRIX(RE_indexes);
  DATA_IVECTOR(nobs_RE);
  DATA_IVECTOR(ln_tau_G_index);
  
  // Derived data 
  int n_RE = nobs_RE.size();

  // Corr matrix 
  DATA_MATRIX(COR);                           // this is the empirical correlation matrix (stacked together)
  
  // INLA features & model configurations (e.g. how to move from mesh to actual data point)      
  DATA_SPARSE_MATRIX(Aobs);                   //Matrix for interpolating points within triangles (for the observation) - used for the mean spatial field
  DATA_SPARSE_MATRIX(Ast);                    //Same but now divided for each year (each observation in a year) - used for the spatio-temporal field
  DATA_IVECTOR(A_spatial_index);              //Vector of stations to match up A_st output
  DATA_FACTOR(year_i);                        //Year index for the stations 
  DATA_INTEGER(Nyear);                        //Number of years (needed if including the spatio-temporal effect)
  DATA_INTEGER(Nmesh);                        //Number of mesh
  DATA_INTEGER(Npred);                        //Number of prediction points
  DATA_INTEGER(do_predict);                   //Doing prediction or not (save time if just exploring) 0 = not, 1 = yes (real scale), 2= yes, link scale
  DATA_INTEGER(corr_str);                     //How to model the correlation structure
  DATA_INTEGER(include_dd);                   //Doing prediction or not (save time if just exploring) 0 = not, 1 = yes
  DATA_INTEGER(calc_se);                      //Calculating SD around P(age) (save time if just exploring) 0 = not, 1 = yes
  DATA_MATRIX(X_proj);                        //The design matrix for the prediction
  DATA_SPARSE_MATRIX(A_proj);                 //The A matrix for the prediction
  DATA_IVECTOR(RE_indexes_proj);              //The random effect index to include in the prediction
  DATA_MATRIX(proj_NA);                       //The location with NA for prediction
  // DATA_VECTOR_INDICATOR(keep, yobs_vec);      // For one-step predictions; https://rdrr.io/cran/TMB/man/oneStepPredict.html
  DATA_ARRAY_INDICATOR(keep, yobs);           // https://rdrr.io/cran/TMB/man/oneStepPredict.html
  DATA_ARRAY(to_keep);                        // https://rdrr.io/cran/TMB/man/oneStepPredict.html
  
  // Distribution
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  DATA_INTEGER(ARorIID);                      // AR=1 or IID=0 for both the presence absence and pos model
  DATA_INTEGER(sim);
  DATA_VECTOR(boxcox_pow);                    // boxcox transformation powcoef component
  
  // Parameters to estimate       
  PARAMETER_MATRIX(beta);                     // coefficient associated with the fiexed effect: Neffect x Nage
  PARAMETER_MATRIX(omega);                    // The mean spatial effects (by age class): Nmesh x Nage (not used for annual barrier model)
  PARAMETER_ARRAY(epsilon_st);                // The spatio-temporal effect: Nmesh x Ntime x Nage
  PARAMETER_VECTOR(transf_rho);               //The autocorrelation value between age or year
  PARAMETER(transf_rho_age);                  //The autocorrelation value between age groups only used when corr_str = 3

  PARAMETER_VECTOR(logKappa);                 //Spatial scale parameter in Matern covariance structures
  PARAMETER_VECTOR(logTauO);                  //Precision parameter for the average spatial field
  PARAMETER_VECTOR(logTauE);                  //Precision parameter for the spatio-temporal variation 
  PARAMETER_VECTOR(logsds);                   //Marginal variance of the species cov matrix

  PARAMETER_VECTOR(thetaf);                   // tweedie only
  PARAMETER_VECTOR(ln_phi);                   // sigma / dispersion / etc.
  
    // Parameters for the random effect (acting on the interecept)
    PARAMETER_MATRIX(ln_tau_G);               // random intercept sigmas
    PARAMETER_MATRIX(RE);                     // random intercept deviations
  
    // Parameters related to the density-dependent effect (logistic model parameters)
    PARAMETER(s50);                             // point of 50% value
    PARAMETER(logslope);                        // steepness of the slope of the logistic curve
    
  // derived parameters
  vector<Type> sds=exp(logsds);
  int N_choice=Nage;
  if(corr_str== 1) N_choice = Nyear;
  if(corr_str == 2) N_choice = Nage;
     
    vector<Type> rho(N_choice);
    for (int t=0; t<N_choice; t++){
      rho(t) =f(transf_rho(t));
    } 

  Type rho_age = f(transf_rho_age);
  vector<Type> range = sqrt(Type(8.0)) / exp(logKappa);
  vector<Type> kappa = exp(logKappa);
  vector<Type> sigma_E = 1 / sqrt(Type(4.0) * M_PI * exp(Type(2.0) * logTauE) * exp(Type(2.0) * logKappa));
  Type slope = exp(logslope);
  
  // ======================== Calculate the linear predictors (link space) then back-transform ========================
  
  matrix<Type> eta(Nobs, Nage);     // this is at the link scale 
  matrix<Type> fixed_noRE(Nobs, Nage);   // this is at the fixed effect
  matrix<Type> fixed(Nobs, Nage);   // this is at the fixed effect
  matrix<Type> mu(Nobs, Nage);      // this is at the real scale 
  
  // Step 1: Add the contribution of the fixed effects
  for(int j=0;j<Nage; j++){
    eta.col(j) = X * vector<Type>(beta.col(j));
    fixed_noRE.col(j) = eta.col(j);
  }

  // Step 2: Now add the contribution of the random effects when existing (on the intercept - no slope):
  if (n_RE > 0){
    for (int i = 0; i < Nobs; i++){
      for (int j = 0; j < Nage; j++){
				int temp = 0;
        for (int k = 0; k < n_RE; k++) {
          if (k == 0) eta(i,j) += RE(RE_indexes(i, k),j); // record it
          if (k > 0) {
            temp += nobs_RE(k - 1);
            eta(i,j) += RE(RE_indexes(i, k) + temp, j); // record it
            fixed_noRE(i,j) += RE(RE_indexes(i, k) + temp, j); // for including RE other than the vessel effect (1st RE variable is always vessel effect)
          }
        }
        fixed(i,j) = eta(i,j);
      }
    } 
  }

  
  if (spatialmodel == 1){  
    // Step 3: Add the contribution of the spatial random field (on the intercept)
    // Here we are "projecting" the spatiotemporal and spatial random effects to the
    // locations of the data using the INLA 'A' matrices, Aobs. Because everything is calculated at the mesh vertice level in INLA
    matrix<Type> omega_A(Nobs, Nage);
    for(int j=0;j<Nage; j++){
      omega_A.col(j) = Aobs* vector<Type>(omega.col(j));
      eta.col(j) = vector<Type>(eta.col(j)) + vector<Type>(omega_A.col(j));
    }
    
    
    // Step 4: Add the contribution of the spatio-temporal effect
    // Same using the A matrices but a little more complicated because interporlation needs to be done for each year's observation separately
    // Begin with calculating the effect at each observation level across year and age
    array<Type> epsilon_st_A(Ast.rows(), Nyear, Nage);
    array<Type> epsilon_st_A_temp(Ast.rows(), Nyear);
    vector<Type> Tmp_st(Nmesh);
    for (int j = 0; j < Nage; j++){
      epsilon_st_A_temp.setZero();
      for (int i = 0; i < Nyear; i++){
        Tmp_st.setZero();
        for (int k=0; k < Nmesh; k++) {
          Tmp_st(k) = epsilon_st(k,i,j);
        }
        epsilon_st_A_temp.col(i) = Ast * Tmp_st;
        for (int k=0; k < Ast.rows(); k++) {
          epsilon_st_A(k,i,j) = epsilon_st_A_temp(k,i);
        }
      }
    }
    // now finding the YEAR the observation takes place and assign that effect value from epsilon_st_A
    matrix<Type> epsilon_st_A_mat(Nobs, Nage);
    epsilon_st_A_mat.setZero();
    for (int i = 0; i < Nobs; i++){
      for (int j = 0; j < Nage; j++){
        epsilon_st_A_mat(i,j) = epsilon_st_A(A_spatial_index(i), year_i(i),j);
        // adding the density dependent effect (if needed)
        if (include_dd == 1) {
          eta(i,j) *= logistic_func(epsilon_st_A_mat(i,j), s50, slope);
        }
        if (include_dd == 0) {
          eta(i,j) += epsilon_st_A_mat(i,j);
        }
      }
    }
    
    REPORT(epsilon_st_A_mat);
    
  }
  
  // Step 5: back transform to the real scale
  for (int i = 0; i < Nobs; i++){
    for (int j = 0; j < Nage; j++){
      Type powcoef = boxcox_pow(j);
      mu(i,j) = InverseLink(eta(i,j), link, powcoef);
    }
  }
  
  
  // ======================== The likelihood components ========================
  
  // Defining the NLL
  Type NLL = 0;
  
  // The random effects (only add the contribution if existing)
  if (n_RE > 0){
    for (int j = 0; j < Nage; j++){
      for (int k = 0; k < n_RE; k++){
        for (int g = 0; g < nobs_RE(k); g++) {
          if (k==0){
            NLL -= dnorm(RE(g,j), Type(0.0), exp(ln_tau_G(ln_tau_G_index(g),j)), true);
            if (sim == 1) SIMULATE{RE(g,j) = rnorm(Type(0), exp(ln_tau_G(ln_tau_G_index(g),j)));}
          }
          if (k>0){
            int temp4 =0;
            for (int bb=0; bb<(k-1); bb++) temp4 += nobs_RE(bb);
            NLL -= dnorm(RE(g+temp4,j), Type(0.0), exp(ln_tau_G(ln_tau_G_index(g+temp4),j)), true);
            if (sim == 1) SIMULATE{RE(g+temp4,j) = rnorm(Type(0), exp(ln_tau_G(ln_tau_G_index(g+temp4),j)));}
          }
        }
      }
    }
  }
  
  // The spatial random effect
  
  if (spatialmodel == 1){  
    
    // The spatio-temporal random effect (same kappa = same spatial range): spatial distribution within a year is correlated between ages 
      Eigen::SparseMatrix<Type> Q;
      if (corr_str == 1){
        array<Type> epsilon_st_temp1(Nmesh, Nage);
        for (int t = 0; t<Nyear; t++){
          Q = Q_spde(spde_barrier, kappa(t), barrier_scaling);
          epsilon_st_temp1.setZero();
          for(int i=0; i<Nage; i++){
            for (int j=0; j<Nmesh; j++){
              epsilon_st_temp1(j,i) = epsilon_st(j,t,i);
            }
          }
          
          if (ARorIID == 1) {
            NLL += SCALE(SEPARABLE(AR1(rho(t)), GMRF(Q,false)), 1.0/exp(logTauE(t)))(epsilon_st_temp1);
            
            if (sim == 1) {
              SIMULATE {SEPARABLE(AR1(rho(t)), GMRF(Q, false)).simulate(epsilon_st_temp1);}
              epsilon_st_temp1 *= 1./exp(logTauE(t));
              for(int i=0; i<Nage; i++){
                for (int j=0; j<Nmesh; j++){
                  epsilon_st(j,t,i) = epsilon_st_temp1(j,i) ;
                }
              }
            }
          }
          
          if (ARorIID == 2) {
            
            array<Type> covar(Nage,Nage,Nyear);
            covar.setZero();
            matrix<Type> cov = COR.block(Nage*t,0,Nage,Nage);
            
            for(int irow=0; irow<Nage; irow++){
              for(int icol=0; icol<Nage; icol++){
                cov(irow,icol) = sds(irow)*sds(icol)*cov(irow,icol);
                covar(irow,icol,t) = cov(irow,icol);
              }
            }
            
            density::MVNORM_t<Type> nldens(cov);
            NLL += SCALE(SEPARABLE(nldens, GMRF(Q,false)), 1.0/exp(logTauE(t)))(epsilon_st_temp1);
            if (sim == 1) {
              SIMULATE {SEPARABLE(nldens, GMRF(Q, false)).simulate(epsilon_st_temp1);}
              epsilon_st_temp1 *= 1./exp(logTauE(t));
              for(int i=0; i<Nage; i++){
                for (int j=0; j<Nmesh; j++){
                  epsilon_st(j,t,i) = epsilon_st_temp1(j,i) ;
                }
              }
            }
            
            // covariance array (Nage x Nage x Nyear)
            REPORT(covar);

          }
          
          if (ARorIID == 0) {
            for(int i=0; i<Nage; i++){
              NLL += VECSCALE(GMRF(Q,false), 1.0/(exp(logTauE(t))*sds))(epsilon_st_temp1.col(i));
            }
            
            if (sim == 1) {
              vector<Type> tmp(epsilon_st_temp1.rows());
              SIMULATE {GMRF(Q, false).simulate(tmp);}
              for(int i=0; i<Nage; i++){
                epsilon_st_temp1.col(i) = tmp/(exp(logTauE(t))*sds(i));
              }
              for(int i=0; i<Nage; i++){
                for (int j=0; j<Nmesh; j++){
                  epsilon_st(j,t,i) = epsilon_st_temp1(j,i) ;
                }
              }
            }
            
          }
          
        }
      }    
        
      // Case when for each age class, spatial distribution is correlated over time
      if (corr_str == 2){
        array<Type> epsilon_st_temp11(Nmesh, Nyear);
        for (int i = 0; i<Nage; i++){
          Q = Q_spde(spde_barrier, kappa(i), barrier_scaling);
          epsilon_st_temp11.setZero();
          for(int t=0; t<Nyear; t++){
            for (int j=0; j<Nmesh; j++){
              epsilon_st_temp11(j,t) = epsilon_st(j,t,i);
            }
          }  
          if (ARorIID == 1) NLL += SCALE(SEPARABLE(AR1(rho(i)), GMRF(Q,false)), 1.0/exp(logTauE(i)))(epsilon_st_temp11);
          if (ARorIID == 0) {
            for(int t=0; t<Nyear; t++){
              NLL += SCALE(GMRF(Q,false), 1.0/exp(logTauE(i)))(epsilon_st_temp11.col(t));
            }
          }
       
          if (sim == 1) {
            vector<Type> tmp(epsilon_st_temp11.rows());
            if (ARorIID == 1) {
              SIMULATE {SEPARABLE(AR1(rho(i)), GMRF(Q, false)).simulate(epsilon_st_temp11);}
              epsilon_st_temp11 *= 1./exp(logTauE(i));
            }
            if (ARorIID == 0) {
              SIMULATE {GMRF(Q, false).simulate(tmp);}
              for(int t=0; t<Nyear; t++){
                epsilon_st_temp11.col(t) = tmp/exp(logTauE(i));
              }
            }  
  					for(int t=0; t<Nyear; t++){
              for (int j=0; j<Nmesh; j++){
                epsilon_st(j,t,i) = epsilon_st_temp11(j,t) ;
              }
            }        
          }
        }
      }
  
      // Case when for spatial distribution is correlated over time and age
      if (corr_str == 3){
        Q = Q_spde(spde_barrier, kappa(0), barrier_scaling);
        NLL += SEPARABLE(AR1(rho_age), SCALE(SEPARABLE(AR1(rho(0)), GMRF(Q,false)), 1.0/exp(logTauE(0))))(epsilon_st);
    
        if (sim == 1) {
          array<Type> epsilon_st_scale(Nmesh, Nyear, Nage);
          SIMULATE {SEPARABLE(AR1(rho_age), SEPARABLE(AR1(rho(0)), GMRF(Q, false))).simulate(epsilon_st_scale);}
          epsilon_st = epsilon_st_scale/exp(logTauE(0));
        }
      }
  }
    
  // The observation likelihood: Only tweedie at the moment
  vector<Type> phi = exp(ln_phi);
  for (int i=0; i<Nobs; i++) {
    for (int j = 0; j < Nage; j++){
      // case tweedie_family:
      Type s1 = invlogit(thetaf(j)) + Type(1.0);
      if (!isNA(yobs(i,j))) NLL -= to_keep(i,j) * keep(i,j) * dtweedie(yobs(i,j), mu(i,j), phi(j), s1, true);
      // if (!isNA(yobs_vec(yobs_index(i,j)))) NLL -= to_keep(i,j) * keep(yobs_index(i,j)) * dtweedie(yobs_vec(yobs_index(i,j)), mu(i,j), phi(j), s1, true);
      if (sim == 1) {
        SIMULATE { 
          yobs(i,j) = rtweedie(mu(i,j), phi(j), s1);
          // yobs_vec(yobs_index(i,j)) = rtweedie(mu(i,j), phi(j), s1);
        }
      }   // break;
    }
  }
  
  
  // ======================== Prediction on new data ========================
  
  if (do_predict == 1) {
    
    array<Type> eta_proj(Npred, Nyear, Nage);   // Combined projection in link space
    array<Type> mu_proj(Npred, Nyear, Nage);    // combined projection in probability scale
    matrix<Type> IA(Nyear, Nage);               // The index of abundance - not corrected (by expert opinion on no-fish zone)
    matrix<Type> IA_corrected(Nyear, Nage);     // The index of abundance corrected for no-fish zone
    matrix<Type> IA_log(Nyear, Nage);           // The index of abundance - log scale
    matrix<Type> IA_corrected_log(Nyear, Nage); // The index of abundance corrected - log-scale
    eta_proj.setZero();
    mu_proj.setZero();
    IA.setZero();
    IA_corrected.setZero();
    
    // the fixed effects
      array<Type> fixed_proj_temp(Npred, Nyear);
      for (int j = 0; j < Nage; j++){
        fixed_proj_temp.setZero();
        for (int i = 0; i < Nyear; i++){
          fixed_proj_temp.col(i) = matrix<Type>(X_proj.block(i*Npred, 0, Npred, beta.rows())) * vector<Type>(beta.col(j));
          for (int k=0; k < Npred; k++) {
            eta_proj(k,i,j) = fixed_proj_temp(k,i);
          }
        }
      }
    
    // The random effects
      
      // Step 2: Now add the contribution of the random effects when existing (on the intercept - no slope):
      if (RE_indexes_proj.size() > 1){
        for (int i = 0; i < Nyear; i++){
          for (int j = 0; j < Nage; j++){
            for (int k=0; k < Npred; k++) {
              int kk = i*Npred+k;
              int temp = 0;
              if (n_RE>1) int temp = nobs_RE(0);
              Type bla = 0;
              if (RE_indexes_proj(kk) != 9999) bla = RE(RE_indexes_proj(kk) + temp, j);
              eta_proj(k,i,j) += bla; // record it
            }
          }
        } 
      }

  array<Type> epsilon_st_A_proj_temp(Npred, Nyear);
  vector<Type> Tmp_st_proj(Nmesh);
  matrix<Type> omega_A_proj(Npred, Nage);
  omega_A_proj.setZero();
   
   if (spatialmodel == 1){  
     
      // the average spatial field
      for(int j=0;j<Nage; j++){
        omega_A_proj.col(j) = A_proj * vector<Type>(omega.col(j));
      }
    
      // the spatio-temporal field
      for (int j = 0; j < Nage; j++){
        epsilon_st_A_proj_temp.setZero();
        for (int i = 0; i < Nyear; i++){
          Tmp_st_proj.setZero();
          for (int k=0; k < Nmesh; k++) {
            Tmp_st_proj(k) = epsilon_st(k,i,j);
          }
          epsilon_st_A_proj_temp.col(i) = A_proj * Tmp_st_proj;
          for (int k=0; k < Npred; k++) {
            eta_proj(k,i,j) += omega_A_proj(k,j) + epsilon_st_A_proj_temp(k,i);
          }
        }
      }
 
   }
       
    // Now back-transform the variable & calculate the index of abundance
      for (int k = 0; k < Nyear; k++){
        for (int j = 0; j < Nage; j++){
          for (int i = 0; i < Npred; i++){
            Type powcoef = boxcox_pow(j);
            mu_proj(i,k,j) = InverseLink(eta_proj(i,k,j), link, powcoef);
            if (!isNA(mu_proj(i,k,j))) IA(k,j) += mu_proj(i,k,j);
            if (!isNA(mu_proj(i,k,j)) & !isNA(proj_NA(i,k))) IA_corrected(k,j) += mu_proj(i,k,j);
          }
          IA_log(k,j) = log(IA(k,j));
          IA_corrected_log(k,j) = log(IA_corrected(k,j));
        }
      }

   
    // Report outputs  
      // REPORT(fixed_proj);               // the predicted fixed effect
      // if (spatialmodel == 1) {
        // REPORT(epsilon_st_A_proj);        // the predicted spatio-temporal changed in spatial field
      // }
      ADREPORT(IA);                     // the predicted IA (not corrected)
      ADREPORT(IA_log);                     // the predicted IA (not corrected)
      ADREPORT(IA_corrected);           // the predicted IA (corrected for no-fish-zone)
      ADREPORT(IA_corrected_log);           // the predicted IA (corrected for no-fish-zone)
      if (calc_se == 0) {
        REPORT(mu_proj);
      }
      if (calc_se == 1) {
        ADREPORT(mu_proj);
      }
    
  }
  
  
  if (do_predict == 2) {
    
    // the fixed effects
    matrix<Type> fixed_proj1(Npred, Nage);
    for (int j = 0; j < Nage; j++){
      fixed_proj1.col(j) = X_proj * vector<Type>(beta.col(j));
    }

    // Report outputs  
    if (do_predict ==2) {
      if (calc_se == 0) {
        REPORT(fixed_proj1);
      }
      if (calc_se == 1) {
        ADREPORT(fixed_proj1);
      }
    }
    
    
  }
  
  
  
  // ======================== Reporting ========================
  
  // Model parameters
  REPORT(beta);
  // REPORT(logTauE);
  // REPORT(sigma_E);
  // REPORT(range);
  
  // Derived parameters on species distribution
  REPORT(fixed_noRE);
  REPORT(fixed);
  // if (spatialmodel == 1){
  //   REPORT(epsilon_st_A_mat); 
  // } 
  // REPORT(epsilon_st_A);

  // REPORT(Q);		
  // REPORT(Q1);
  // REPORT(Q3);		
  if (sim == 1) {
    SIMULATE {
    REPORT(yobs); 
    // REPORT(yobs_vec); 
    REPORT(RE);
    if (spatialmodel == 1) REPORT(epsilon_st);
    }
  }
  // Some outputs for the diagnostics
  REPORT(mu);
  
  
  return NLL;
  
  
}
