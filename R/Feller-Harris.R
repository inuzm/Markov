# Simulate from Markov processes with transition functions given by
# P_t(x, A) = (1 - g(t)) Q(A) + g(t) delta_x(A)
# with Q a probability measure and delta_x the Dirac delta.

# First algorithm. It can be seen that this process can be obtained via a
# Poisson process.

# n is the number of steps
# t is a vector of increasing times of length n+1 (starting with 0)
# lambda_a is the rate for the Poisson process.
# theta_a is a weigth…
# x0 is the starting point (possibly not given)
# dist is the distribution Q
# ... parameters for dist
algy <- function(n, t, lambda_a = 1, theta_a = 1, x0 = NULL, dist, ...){
    X <- numeric(n+1)
    X[1] <- ifelse(is.null(x0), dist(n = 1, ...), x0)
    for(i in 1:n){
        aux1 <- rpois(n = 1, lambda = lambda_a * (t[i+1] - t[i]))
        aux2 <- rbinom(n = 1, size = 1, prob = (1+theta_a)^(-aux1))
        X[i+1] <- ifelse(
            aux2 == 1,
            X[i],
            dist(n = 1, ...)
        )
    }
    return(X)
}

# Second algorithm. Directly from the the transition function that satisfies
# Chapman–Kolmogorov equations.

# n is the number of steps
# t is a vector of increasing times of length n+1 (starting with 0)
# alpha_a is the parameter given in the transition
# x0 is the starting point (possibly not given)
# dist is the distribution Q
# ... parameters for dist
algz <- function(n, t, alpha_a = 1, x0 = NULL, dist, ...){
    X <- numeric(n+1)
    X[1] <- ifelse(is.null(x0), dist(n = 1, ...), x0)
    for(i in 1:n){
        aux <- rbinom(n = 1, size = 1, prob = exp(-alpha_a * (t[i+1] - t[i])))
        X[i+1] <- ifelse(
            aux == 1,
            X[i],
            dist(n = 1, ...)
        )
    }
    return(X)
}
