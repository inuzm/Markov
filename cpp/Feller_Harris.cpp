#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
// First algorithm. It can be seen that this process can be obtained via a
// Poisson process
NumericMatrix feller_harris(NumericVector tspan, double lambda_a, double theta_a, NumericVector x0, int nsim, Rcpp::Function qfun) {
  //
  // INPUT
  // nsim .- Number of different simulations of process
  // t    .- Vector of length n with times of simulations with t(0) = 0;
  // lambda_a is the rate for the Poisson process.
  // theta_a is a weigth…
  // x0 is the starting point
  //
  // OUTPUT
  // X .- Matrix (dimensions: nsim x n) with simulated entries along tspan of the process
  //      where n = length(tspan).
  //
  // WARNINGS
  // Can get memory error if length(t)*nsim is too large.
  // For large values of theta_a, lambda_a, might throw NaN
  //
  
  //Get number of points to run the simulation
  int n = tspan.size();
  
  //Instantiate vectors
  NumericVector aux(nsim);
  NumericVector ones(nsim, 1.0);
  NumericVector Q(nsim);
  double dt; 
  
  //Create simulation matrix
  NumericMatrix X(nsim, n);   
  
  //First entry as x0
  X(_, 0) = x0;
  
  for (int i = 0; i < (n-1); i++) {
      
      //Create auxiliary variables
      dt       = tspan(i + 1) - tspan(i);
      for (int j = 0; j < nsim; j++){
        aux(j) = R::rbinom(1, pow(1.0 + theta_a, -R::rpois(lambda_a*dt)));
      }
      
      Q = qfun();
      X(_, i + 1) = aux*X(_, i) + (ones - aux)*Q;
  }
  
  return X;
  
}

// [[Rcpp::export]]
// Second algorithm. Directly from the the transition function that satisfies
// Chapman–Kolmogorov equations.
NumericMatrix feller_harris_ck(NumericVector tspan, double alpha_a, NumericVector x0, int nsim, Rcpp::Function qfun) {
  //
  // INPUT
  // nsim    .- Number of different simulations of process
  // t       .- Vector of length n with times of simulations with t(0) = 0;
  // alpha_a .- the parameter given in the transition
  // x0      .- is the starting point
  //
  // OUTPUT
  // X .- Matrix (dimensions: nsim x n) with simulated entries along tspan of the process
  //      where n = length(tspan).
  //
  // WARNINGS
  // Can get memory error if length(t)*nsim is too large.
  // For large values of alpha, might throw NaN
  //
  
  //Get number of points to run the simulation
  int n = tspan.size();
  
  //Instantiate vectors
  NumericVector aux(nsim);
  NumericVector ones(nsim, 1.0);
  NumericVector Q(nsim);
  double dt; 
  
  //Create simulation matrix
  NumericMatrix X(nsim, n);   
  
  //First entry as x0
  X(_, 0) = x0;
  
  for (int i = 0; i < (n-1); i++) {
    
    //Create auxiliary variables
    dt       = tspan(i + 1) - tspan(i);
    aux      = rbinom(nsim, 1, exp(-alpha_a * dt));
    Q        = qfun();
    X(_, i + 1) = aux*X(_, i) + (ones - aux)*Q;
  }
  
  return X;
  
}

/*** R
rfeller.harris <- function(tspan, lambda_a = 1, theta_a = 1, alpha_a = lambda_a*theta_a/(1 + theta_a), 
                           nsim = 1, x0 = rnorm(nsim, mean = 0, sd = 1),
                           type = c("poissonprocess","chapmankolmogorov"), dist, ...){
  
#Controls for simulation
  if (!is.vector(tspan)){ 
    stop("t must be a time vector") 
  }
  if (tspan[1] != 0){ 
    stop("t must start at time t[1] = 0") 
  }
  if (lambda_a < 0){ 
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
  
  #Check if distribution is specified; if not, set to rnorm
  if (hasArg(dist)){
    inputdistribution <- function(){dist(n = nsim, ...)}  
  } else {
    inputdistribution <- function(){rnorm(nsim)}
  }
  
  simulation <- switch(type[1],
                       poissonprocess    = feller_harris(tspan, lambda_a, theta_a, x0, nsim, inputdistribution),
                       chapmankolmogorov = feller_harris_ck(tspan, alpha_a, x0, nsim, inputdistribution),
                       stop("Invalid simulation type"))
    
    
#Run simulation
    return(simulation)
}
*/
