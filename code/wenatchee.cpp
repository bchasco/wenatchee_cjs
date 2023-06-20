#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data
  DATA_ARRAY(c_it);     // c_it history matrix
  DATA_ARRAY(z_it);     // c_it history matrix
  DATA_IVECTOR(f_i); // Number of survey locations
  DATA_VECTOR(t_i);        // Time when fish leave Wenatchee mainstem.
  // DATA_INTEGER(n_time);     // Number of time points
  
  // Parameters
  PARAMETER_VECTOR(f_phi); // Logit-transformed c_it probabilities
  PARAMETER_VECTOR(f_p); // Logit-transformed detection probabilities
  PARAMETER(lam);     // Log lambda for time-density model
  
  int n_ind = c_it.dim[0]; //numver of individuals 
  int n_loc = c_it.dim[1]; //number of locations

  // Initialize log-likelihood
  Type nll = 0;
  

  // Loop over individuals
  for (int i = 0; i < n_ind; i++) {
    
    //First observation
    Type eta_phi_i = f_phi[0] + lam * t_i(i);
    Type phi_i = exp(eta_phi_i)/(1+exp(eta_phi_i)); // c_it probabilities

    Type eta_p_i = f_p(0);
    Type p_i = exp(eta_p_i)/(1+exp(eta_p_i));
    
    //likelihood of the first observation
    nll -= log(dbinom(Type(1.),Type(1.),phi_i));
    //Observed process  
    // nll -= log(dbinom(c_it(i,f_i[i]-1),z_it(i,j),p_i));
    
    // Loop over survey locations
    for (int j = f_i[i]; j < n_loc; j++) {
      // Loop over capture histories 

      //Detection probability
      // Type eta_p_i = f_p(0);
      // Type p_i = exp(eta_p_i)/(1+exp(eta_p_i));

      //Biological process
      nll -= log(dbinom(z_it(i,j-1),Type(1.),phi_i));
      
      //Observed process  
      nll -= log(dbinom(z_it(i,j),c_it(i,j),p_i));
      
        // // Calculate the log-likelihood contribution for this cap occasion
        // if (c_it(i, j) == 1)
        //   nll -= log(p_capt); // Individual was capd
        // else
        //   nll -= log(1 - p_capt); // Individual was not capd
        // 
        // // Calculate the log-likelihood contribution for this detection occasion
        // if (j > 0 && c_it(i, j - 1) == 1)
        //   nll -= log(p_detect); // Individual was detected
        // else
        //   nll -= log(1 - p_detect); // Individual was not detected
    }
  }
  
  return nll;
}

