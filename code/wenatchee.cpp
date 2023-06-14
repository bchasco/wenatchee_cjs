#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX(capture);     // Capture history matrix
  DATA_INTEGER(n_locations); // Number of survey locations
  DATA_VECTOR(time);        // Time points
  DATA_INTEGER(n_time);     // Number of time points
  
  // Parameters
  PARAMETER_VECTOR(logit_p); // Logit-transformed capture probabilities
  PARAMETER_VECTOR(logit_d); // Logit-transformed detection probabilities
  PARAMETER(log_lambda);     // Log lambda for time-density model
  
  // Transform parameters
  vector<Type> p = exp(logit_p) / (1 + exp(logit_p)); // Capture probabilities
  vector<Type> d = exp(logit_d) / (1 + exp(logit_d)); // Detection probabilities
  
  // Initialize log-likelihood
  Type nll = 0;
  
  // Loop over individuals
  for (int i = 0; i < capture.rows(); i++) {
    int n_captures = capture.cols(); // Number of capture occasions for the individual
    
    // Loop over survey locations
    for (int j = 0; j < n_locations; j++) {
      // Loop over capture occasions
      for (int k = 0; k < n_captures; k++) {
        // Probability of capture at location j and time k
        Type logit_p_capt = logit_p[j] + log_lambda * time[k];
        Type p_capt = exp(logit_p_capt) / (1 + exp(logit_p_capt));
        
        // Probability of detection at location j
        Type p_detect = d[j];
        
        // Calculate the log-likelihood contribution for this capture occasion
        if (capture(i, k) == 1)
          nll -= log(p_capt); // Individual was captured
        else
          nll -= log(1 - p_capt); // Individual was not captured
        
        // Calculate the log-likelihood contribution for this detection occasion
        if (k > 0 && capture(i, k - 1) == 1)
          nll -= log(p_detect); // Individual was detected
        else
          nll -= log(1 - p_detect); // Individual was not detected
      }
    }
  }
  
  return nll;
}

