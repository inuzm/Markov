#Plot simulations
plotsims <- function(tspan, simmat, xaxis = TeX("Time ($t_k$)"), yaxis = TeX("Simulations of $X_{t_k}$"), 
                     mytitle = "A Stochastic Process", mycolors = rainbow(nrow(simmat))){
  
  #Each row in simmat is one simulation along tspan
  if (length(tspan) != ncol(simmat)){
    stop("Error: Time span, tspan, does not correspond to matrix simulation period.")
  }
  
  #Calculate number of rows and columns
  nsim <- nrow(simmat)
  n    <- length(tspan)
  
  #Obtain min and max for plot
  pmin <- Inf
  pmax <- -Inf
  
  #Create empty ggplot object
  simplot <- ggplot()
  
  for (i in 1:nsim){
    pmin <- min(min(simmat[i,]), pmin)
    pmax <- max(max(simmat[i,]), pmax)
    simplot <- simplot + geom_line(aes(x = x, y = y), 
                                   data = data.frame(x = tspan, y = simmat[i,]), color = mycolors[i])
  }
  
  #Add labels
  simplot <- simplot + xlab(xaxis) + ylab(yaxis) + theme_classic() + ylim(pmin, pmax)
  
  #Add empirical stationary distribution plot
  dplot <- ggplot() + 
    geom_density(aes(x = x), kernel = "gaussian", color = "deepskyblue1",
                 data = data.frame(x = simmat[,ncol(simmat)])) + xlim(pmin, pmax) +
    theme_classic() + xlab(TeX('Distribution of $X_{t_n}$')) + ylab("Kernel density approximation") +  
    coord_flip()
  
  #Arrange in grid
  mygrid <- grid.arrange(simplot, dplot, nrow = 1, widths = c(2,1), heights = c(1), top = mytitle)
  
  
  return(mygrid)
  
}