# Simulate from an Ornstein–Uhlenbeck Process seeing it as a Markov Process.

# First algorithm. Direct from transition functions.
ornstein.uhlenbeck <- function(n, t, mu = 0, tau = 1, alpha = 1){
    X <- numeric(n+1)
    X[1] <- rnorm(n = 1, mean = mu, sd = sqrt(tau))
    for(i in 1:n){
        X[i+1] <- rnorm(
            n = 1,
            mean = X[i] * exp(-alpha * t) + mu * (1 - exp(-alpha * t)),
            sd = sqrt(tau * (1 - exp(-2 * alpha * t)))
        )
    }
    return(X)
}

# Second algorithm, seeing that the transition function can be seen as a mixture.
ornstein.uhlenbeck.aum <- function(n, t, mu = 0, tau = 1, alpha = 1){
    X <- numeric(n+1)
    X[1] <- rnorm(n = 1, mean = mu, sd = sqrt(tau))
    for(i in 1:n){
        aux <- rnorm(
            n = 1,
            mean = X[i],
            sd = sqrt(tau * (exp(alpha * t) - 1))
        )
        X[i+1] <- rnorm(
            n = 1,
            mean = aux * exp(-alpha * t) + mu * (1 - exp(-alpha * t)),
            sd = sqrt(tau * (1 - exp(- alpha * t)))
        )
    }
    return(X)
}

# Third algorithm, discretizing the stochastic differential equation associated with
# the transition function given in the first algorithm
ornstein.uhlenbeck.ede <- function(n, t, mu = 0, tau = 1, alpha = 1, x0 = rnorm(n = 1, mean = mu, sd = sqrt(tau))){
    X <- numeric(n+1)
    X[1] <- x0
    for(i in 1:n){
        aux <- rnorm(
            n = 1,
            mean = 0,
            sd = sqrt(t[i+1] - t[i])
        )
        X[i+1] <- X[i] + alpha * (mu - X[i]) * (t[i+1] - t[i]) + sqrt(2 * alpha * tau) * aux
    }
    return(X)
}
