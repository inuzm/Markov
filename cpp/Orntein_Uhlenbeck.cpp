// NOTES
// I always get lost with random number generation in Rcpp. I'm not including
// an option to set seed for now. This can be done as established here:
// http://gallery.rcpp.org/articles/random-number-generation/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// First algorithm. Direct from transition functions.
NumericMatrix ornstein_uhlenbeck(NumericVector tspan, double mu, double tau, double alpha, int nsim, NumericVector x0) {
  //
  // INPUT
  // nsim .- Number of different simulations of process
  // t    .- Vector of length n with times of simulations with t(0) = 0;
  // mu, tau, alpha .- Parameters of process
  //
  // OUTPUT
  // X .- Matrix (dimensions: nsim x n) with simulated entries along tspan of the process
  //      where n = length(tspan).
  //
  // WARNINGS
  // Can get memory error if length(t)*nsim is too large.
  // For large values of tau, alpha, might throw NaN
  //
  
  //Get number of points to run the simulation
  int n = tspan.size();
  
  //Create simulation matrix
  NumericMatrix X(nsim, n); 
  
  //Create numeric vector
  NumericVector Z(nsim);
  NumericVector mean(nsim); 
  double desvest; 
  double dt; 
  
  //For first entry simulate nsim random normals
  X(_, 0) = x0;
  for (int i = 0; i < (n-1); i++) {
    
     //Create auxiliary variables
     dt       = tspan(i + 1) - tspan(i);
     Z        = rnorm(nsim, 0, 1);                                       //Standard normal
     desvest  = sqrt( tau * (1 - exp(-2 * alpha * dt)) );                //Variance of transition 
     mean     = X(_,i) * exp(-alpha * dt) + mu * (1 - exp(-alpha * dt)); //Mean of transition
     
     //Simulate normal transition
     X(_, i + 1) = mean + desvest*Z;
     
  }
  
  return X;
}

// [[Rcpp::export]]
// Second algorithm, seeing that the transition function can be seen as a mixture.
NumericMatrix ornstein_uhlenbeck_aum(NumericVector tspan, double mu, double tau, double alpha, int nsim, NumericVector x0) {
  //
  // INPUT
  // nsim .- Number of different simulations of process
  // t    .- Vector of length n with times of simulations with t(0) = 0;
  // mu, tau, alpha .- Parameters of process
  //
  // OUTPUT
  // X .- Matrix (dimensions: nsim x n) with simulated entries along tspan of the process
  //      where n = length(tspan).
  //
  // WARNINGS
  // Can get memory error if length(t)*nsim is too large.
  // For large values of tau, alpha, might throw NaN
  //
  
  //Get number of points to run the simulation
  int n = tspan.size();
  
  //Create simulation matrix
  NumericMatrix X(nsim, n); 
  
  //Create numeric vector
  NumericVector Z1(nsim);
  NumericVector Z2(nsim);
  NumericVector mean(nsim); 
  NumericVector aux(nsim);
  double desvest; 
  double dt; 
  
  //For first entry simulate nsim random normals
  X(_, 0) = x0;
  for (int i = 0; i < (n-1); i++) {
    
    //Create auxiliary variables
    dt       = tspan(i + 1) - tspan(i);
    Z1       = rnorm(nsim, 0, 1);      //Standard normal for aux
    Z2       = rnorm(nsim, 0, 1);      //Standard normal for X[i]
    aux      = X(_,i) + desvest*Z1;
    desvest  = sqrt( tau * (exp(alpha * dt) - 1) );                  //Variance of transition 
    mean     = aux * exp(-alpha * dt) + mu * (1 - exp(-alpha * dt)); //Mean of transition
    
    //Simulate normal transition
    X(_, i + 1) = mean + desvest*Z2;
    
  }
  
  return X;
}

// [[Rcpp::export]]
// Third algorithm, discretizing the stochastic differential equation associated with
// the transition function given in the first algorithm
NumericMatrix ornstein_uhlenbeck_ede(NumericVector tspan, double mu, double tau, double alpha, int nsim, NumericVector x0) {
  //
  // INPUT
  // nsim .- Number of different simulations of process
  // t    .- Vector of length n with times of simulations with t(0) = 0;
  // x0   .- Initial vector distribution
  // mu, tau, alpha .- Parameters of process
  //
  // OUTPUT
  // X .- Matrix (dimensions: nsim x n) with simulated entries along tspan of the process
  //      where n = length(tspan).
  //
  // WARNINGS
  // Can get memory error if length(t)*nsim is too large.
  // For large values of tau, alpha, might throw NaN
  // Use only where abs(tspan(i) - tspan(i+1)) <= 0.1
  
  //Get number of points to run the simulation
  int n = tspan.size();
  
  //Create simulation matrix
  NumericMatrix X(nsim, n); 
  
  //Create numeric vector
  NumericVector Z(nsim);
  double dt; 
  
  //For first entry simulate nsim random normals
  X(_, 0) = x0;
  for (int i = 0; i < (n-1); i++) {
    
    //Create auxiliary variables
    dt       = tspan(i + 1) - tspan(i);
    Z        = rnorm(nsim, 0, 1);        
    
    //Simulate normal transition
    X(_, i + 1) = X(_, i) + alpha*(mu - X(_, i))*dt + sqrt(2 * alpha * tau * dt)*Z;
    
  }
  
  return X;
}

/*** R
rornstein.uhlenbeck <- function(tspan, mu = 0, tau = 1, alpha = 1, nsim = 1, type = c("direct","mixture","sde"),
                                x0 = rnorm(nsim, mean = mu, sd = sqrt(tau))){
  
  #Controls for simulation
  if (!is.vector(tspan)){ 
    stop("t must be a time vector") 
  }
  if (tspan[1] != 0){ 
    stop("t must start at time t[1] = 0") 
  }
  if (tau < 0){ 
    stop("tau must be positive") 
  }
  if (nsim < 1){ 
    stop("nsim is number of simulations must be > 1") 
  }
  if (tspan[2] - tspan[1] > 0.1){ 
    stop("Discretization only works with small dt.") 
  }
  
  #Ensure n is integer
  nsim <- ceiling(nsim)
  
  simulation <- switch(type[1],
                       direct  = ornstein_uhlenbeck(tspan, mu, tau, alpha, nsim, x0),
                       mixture = ornstein_uhlenbeck_aum(tspan, mu, tau, alpha, nsim, x0),
                       sde     = ornstein_uhlenbeck_ede(tspan, mu, tau, alpha, nsim, x0),
                       stop("Invalid simulation type"))
  
  #Run simulation
  return(simulation)
}



*/
