
# load libraries
library(ggplot2) #library for plotting
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

#load concentration data
df = tibble(
  time = 0:100,
  
)
names(df)=c("time","ca","cb","cc")

# rate function
rxnrate=function(t,c,parms){
  
  # rate constant passed through a list called parms
  k1=parms$k1
  k2=parms$k2
  
  # c is the concentration of species
  
  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1]=-k1*c["A"] #dcA/dt
  r[2]=k1*c["A"]-k2*c["B"] #dcB/dt
  r[3]=k2*c["B"] #dcC/dt
  
  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))
  
}

# predicted concentration for a given parameter set
cinit=c(A=1,B=0,C=0)
parms=list(k1=2,k2=1)
t = seq(0, 5, by=0.1)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
plot(out)
df = as_tibble(as.matrix(out)) %>% 
  mutate(
    time = as.numeric(time),
    A = as.numeric(A) + rnorm(nrow(.), 0, 0.05),
    B = as.numeric(B) + rnorm(nrow(.), 0, 0.05),
    C = as.numeric(C) + rnorm(nrow(.), 0, 0.05),
    )

ssq=function(parms){
  
  # inital concentrations
  cinit=c(A=1,B=0,C=0)
  # time points for which conc is reported
  # include the points where data is available
  t=c(seq(0,5,0.1),df$time)
  t=sort(unique(t))
  # parameters from the parameter estimation routine
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% df$time,]
  # Evaluate predicted vs experimental residual
  preddf= outdf %>% 
    pivot_longer(-time, values_to = "conc", names_to = "species") %>% 
    filter(time %in% df$time)
  
  expdf = df %>% 
    pivot_longer(-time, values_to = "conc", names_to = "species")
  
  ssqres=preddf$conc-expdf$conc
  
  # return predicted vs experimental residual
  return(ssqres)
  
}

# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
parms=list(k1=0.5,k2=0.5)
# fitting
fitval=nls.lm(par=parms,fn=ssq)

# fitting
summary(fitval)
