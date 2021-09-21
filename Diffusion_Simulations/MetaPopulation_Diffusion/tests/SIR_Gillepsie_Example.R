# Simple SIR Simulation where the Gillepsie Algorythm Is Implemented
# Jonathan H. Morgan Using Resources Provided by Ben Bolker (2003, 2015)
# 20 July 2021

# LambertW Implmentation: https://stat.ethz.ch/pipermail/r-help/2003-November/042793.html
# Basic SIR Model with the Gillespie Algorithm Implemented: https://rpubs.com/bbolker/SIRgillespie

# Clear Out Console Script
  cat("\014")

# Setting Work Directory
  setwd("~/Dropbox/My Mac (Jonathan’s MacBook Pro)/Desktop/DNAC/IDEANet/Data_Scripts")
  getwd()

# Options
  options(stringsAsFactors = FALSE)

#################
#   FUNCTIONS   #
#################
  
  # LambertW Arguments
  # {z} (complex) vector of values for which to compute the function
  # {b} (integer) vector of branches: b=0 specifies the principal branch, 0 and -1 are the ones that can take non-complex values
  # {maxiter} maximum numbers of iterations for convergence
  # {eps} convergence tolerance
  # {min.imag} maximum magnitude of imaginary part to chop when returning solutions

  # Nici Schraudolph's lambertw implementation
    lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps, min.imag=1e-9) {
      if (any(round(Re(b)) != b))
        stop("branch number for W must be an integer")
      if (!is.complex(z) && any(z<0)) z=as.complex(z)
      
      w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
 
      v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
      v = v - log(v + as.numeric(v==0))
 
      c = abs(z + exp(-1));
      c = (c > 1.45 - 1.1*abs(b));
      c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
      w = (1 - c)*w + c*v
    
      # Halley iteration
        for (n in 1:maxiter) {
          p = exp(w)
          t = w*p - z
          f = (w != -1)
          t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
          w = w - t
          if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w))) && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
            break
        }
      
      if (n==maxiter) warning(paste("iteration limit (",maxiter,") reached, result of W may be inaccurate",sep=""))
      if (all(Im(w)<min.imag)) w = as.numeric(w)
      return(w)
    }

  # Examples of the Function
    curve(lambertW(x),from=0,to=10, bty='n', col="brown", family='HersheySerif', las=1)
    grid(lwd = 2)
    
    pvec = seq(-1,1,length=40)
    m = outer(pvec, pvec, function(x,y) Re(lambertW(x+y*1i)))
    persp(pvec,pvec,m, theta=290,shade=0.5,zlab="lambertW")
    
    num1 = uniroot(function(x) {x*exp(x)-1},lower=0,upper=1,tol=1e-9)
    abs(lambertW(1)-num1$root)<1e-9
    
  # Functions for computing the event rates and the transitions to be executed when the events occur
    
    # Function Specifying the Rates at Infection and Recovery Occurr
      ratefun <- function(x, p) {
        with(as.list(c(x, p)),{
          c(inf = beta*S*I/N,  ## scale inf by pop size
          recover = gamma*I)
        })
      }
  
    # Function Specifying the State Transitions Based on Counts of Susceptible and Infected
      transfun <- function(x,w) {
        switch(w,
               x + c(-1,1),   ## infection: S-1, I+1
               x + c(0,-1))   ## removal: I-1
      }     
    
  # A Wrapper Function to Run the Simulation with Specified Parameters/Starting Values & Return 
  # Either the Ending State or a Matrix of Event Times and Types
    run <- function(p=c(beta=2, gamma=1, N=100), I0=1, itmax=1e5, ret=c("final","all")) {
                    ret <- match.arg(ret)
                    if (ret=="all") {
                      rmat <- matrix(NA,nrow=itmax,ncol=2,
                      dimnames=list(NULL,c("t","trans")))
                    }
                    x <- c(S=unname(p["N"])-I0, I=I0)       # Number of Susceptible and Infected (99 and 1) at time 1
                    it <- 1                                 # Iteration value at time 1
                    t <- 0                                  # Count of Recovered
                    trans <- c(0,0)                         # Transitions Movement from Susceptible (S-1) to Infected (I + 1)
                    while (x["I"] > 0 & it < itmax) {
                      r <- ratefun(x,p)                     # Rates for Susceptible and Infected based on beta, gamma, N values
                      t <- t + rexp(1, rate=sum(r))         # t = t + random draw of one from an exponential distribution with a rate of sum(r)
                      w <- sample(length(r), size=1,prob=r) # Sample one of 2 values (susceptible or infected events) based on the probabilities of the two.
                      x <- transfun(x,w)                    # Update the number of Suceptible and Infected based on the event that occurred in the last step.
                      if (ret=="all") 
                        rmat[it,] <- c(t,w)
                        it <- it+1
                    }
                    if (ret=="all") return(rmat[!is.na(rmat[,1]),])
                    return(c(S=unname(x["S"]),t=t,it=it))
    }
    
  # Analytic Computation of Expected Final Size based on the Lambert W Function
    finalsize <- function(R0) {
      1 + 1/R0*lambertW(-R0*exp(-R0))
    }
  
###################
#   SIMULATIONS   #
###################
    
  # Simulation: Basic Examples
  # Two Infections, Three Recoveries
    set.seed(101)
    ex0 <- run(ret="all")
    plot(0, type='n', xlab='t', ylab='Transitions', xlim=c(0, 1), ylim=c(1, 2), cex.axis=1.3, family='HersheySerif', las=1, main='', bty='n')
    grid(lwd = 2)
    points(x=ex0[,1], y=ex0[,2], col="brown", pch=16, cex=1.3)
    title('Simulation 1', family='serif', cex.main=2, line=1.25)
    
  # One Infection, One Recovery
    set.seed(101)
    ex0 <- run(p=c(beta=1.1, gamma=1,N=100), ret="all")
    plot(0, type='n', xlab='t', ylab='Transitions', xlim=c(0, 1), ylim=c(1, 2), cex.axis=1.3, family='HersheySerif', las=1, main='', bty='n')
    grid(lwd = 2)
    points(x=ex0[,1], y=ex0[,2], col="brown", pch=16, cex=1.3)
    title('Simulation 1', family='serif', cex.main=2, line=1.25)
    abline(v=0.1761341, col='brown', lty=2)
    
  # Exercises 1 and 2
    
    # Beta is at the Default Value of 2
      trials <- vector('list', 1000)
      for (i in seq_along(trials)) {
        trials[[i]] <- as.numeric(run())
      }
      ex1 <- do.call("rbind", trials)
      ex1 <- as.data.frame(cbind(seq(1, 1000, 1), ex1))
      colnames(ex1) <- c('n', 'S', "t", 'it')
    
    # Beta is specified to be 1.1
      for (i in seq_along(trials)) {
        trials[[i]] <- as.numeric(run(p=c(beta=1.1,gamma=1,N=100)))
      }
      ex2 <- do.call("rbind", trials)
      ex2 <- as.data.frame(cbind(seq(1, 1000, 1), ex2))
      colnames(ex2) <- c('n', 'S', "t", 'it')
      rm(trials)
      
######################
#   VISUALIZATIONS   #
######################
    
    # Plotting Simulation 1 and Simulation 2 Results
      layout.matrix <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
      layout(mat = layout.matrix,
             heights = c(2, 2), # Heights of the two rows
             widths = c(2, 2)) # Widths of the two columns
      layout.show(4)
      
      par(mar = c(4,6,3,3),  family='HersheySerif')
      hist(ex1$S, breaks = 35, main=" ", ylim=c(0, 400), xlab='', ylab='', las=1, cex.axis=1.3, family='HersheySerif')
      mtext(side = 1, text = 's', col = "black", line = 3, cex = 1.5, family='HersheySerif')
      mtext(side = 2, text = "Count", col = "black", line = 4, cex = 1.5, family='HersheySerif')
      title('Simulation 1: S', family='serif', cex.main=1.5, line=1.25)
      
      hist(ex1$t, breaks = 35, main=" ", xlab=' ', ylab=' ', ylim=c(0, 400), las=1, cex.axis=1.3, family='HersheySerif')
      mtext(side = 1, text = 't', col = "black", line = 3, cex = 1.5, family='HersheySerif')
      mtext(side = 2, text = "Count", col = "black", line = 4, cex = 1.5, family='HersheySerif')
      title('Simulation 1: t', family='serif', cex.main=1.5, line=1.25)
      
      hist(ex2$S, breaks = 35, main=" ", ylim=c(0, 500), xlab='', ylab='', las=1, cex.axis=1.3, family='HersheySerif')
      mtext(side = 1, text = 's', col = "black", line = 3, cex = 1.5, family='HersheySerif')
      mtext(side = 2, text = "Count", col = "black", line = 4, cex = 1.5, family='HersheySerif')
      title('Simulation 2: S', family='serif', cex.main=1.5, line=1.25)
      
      hist(ex2$t, breaks = 35, main=" ", xlab='', ylab=' ', ylim=c(0, 500), las=1, cex.axis=1.3, family='HersheySerif')
      mtext(side = 1, text = 't', col = "black", line = 3, cex = 1.5, family='HersheySerif')
      mtext(side = 2, text = "Count", col = "black", line = 4, cex = 1.5, family='HersheySerif')
      title('Simulation 2: t', family='serif', cex.main=1.5, line=1.25)
      
  # Estimating Expected Sizes
    finalsize(2)
    finalsize(1.1)
    