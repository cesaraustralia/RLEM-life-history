
# load libraries
library(ggplot2) #library for plotting
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm

#load concentration data
df = tibble(
  time = 0:100,
  
)
# rate function
rxnrate=function(t,state,parms){
  
    with(as.list(c(state, parms)), {
      # dEd = 0
      # dEw = eggrate * M5 / 2 - t1 * Ew
      Ew = get_Ew(t)
      # Ew = 100
      dM1 = t1 * Ew - t2 * M1 - m1*M1
      dM2 = t2 * M1 - t3 * M2 - m2*M2 
      dM3 = t3 * M2 - t4 * M3 - m3*M3
      dM4 = t4 * M3 - t5 * M4 - m4*M4
      dM5 = t5 * M4 - m5*M5
      list(c(dM1, dM2, dM3, dM4, dM5))
    })
  
}

# predicted concentration for a given parameter set
cinit= c(M1 = 0, M2 = 0, M3 = 0, M4 = 0, M5 = 0)
parms=c(
  t1=0.15, 
  t2=0.15, 
  t3=0.15, 
  t4=0.15, 
  t5=0.15,
  m1=0.05, 
  m2=0.05, 
  m3=0.05, 
  m4=0.05, 
  m5=0.2)

t = seq(100, 300, by=1)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)

df = as_tibble(as.matrix(out)) %>% 
  mutate(
    time = as.numeric(time),
    M1 = as.numeric(M1) + rnorm(nrow(.), 0, 2),
    M2 = as.numeric(M2) + rnorm(nrow(.), 0, 2),
    M3 = as.numeric(M3) + rnorm(nrow(.), 0, 2),
    M4 = as.numeric(M4) + rnorm(nrow(.), 0, 2),
    M5 = as.numeric(M5) + rnorm(nrow(.), 0, 2),
    )

dfpred = df %>%  
  pivot_longer(-time)

dfpred %>% 
  ggplot() + 
  geom_line(aes(time, value, color=name))


##### real data ##### 
id = 1
df = read_csv('./data/clean_data.csv') %>% 
  filter(group_id==id) %>% 
  select(time, M1, M2, M3, M4, M5)

Ew = read_csv('./data/clean_data.csv') %>% 
  filter(group_id==id) %>% 
  select(time, Ew)

get_Ew = function(t) {
  Ew_ind = which(abs(Ew$time - t) == min(abs(Ew$time - t)))
  mean(Ew$Ew[Ew_ind], na.rm=TRUE)
}


df %>% 
  left_join(Ew) %>% 
  pivot_longer(-time) %>% 
  ggplot() + 
  geom_point(aes(time, value, color=name)) +
  geom_line(data=dfpred, aes(time, value, color=name)) 
dfpred
get_Ew(200)

rxnrate=function(t,state,parms){
  
  
  with(as.list(c(state, parms)), {
    # dEd = 0
    # dEw = eggrate * M5 / 2 - t1 * Ew
    Ew = get_Ew(t)
    # Ew = 100
    dM1 = t1 * Ew - t2 * M1 - m1*M1
    dM2 = t2 * M1 - t3 * M2 - m2*M2 
    dM3 = t3 * M2 - t4 * M3 - m3*M3
    dM4 = t4 * M3 - t5 * M4 - m4*M4
    dM5 = t5 * M4 - m5*M5
    list(c(dM1, dM2, dM3, dM4, dM5))
  })
  
}

d0 = read_csv('./data/clean_data.csv')

ssq=function(parms){
  
  SSQRES = c() 
  
  for(i in 1:50){
    df = d0 %>% 
      filter(group_id==id) %>% 
      select(time, M1, M2, M3, M4, M5)
    
    Ew = d0 %>% 
      filter(group_id==id) %>% 
      select(time, Ew)
    
    get_Ew = function(t) {
      Ew_ind = which(abs(Ew$time - t) == min(abs(Ew$time - t)))
      mean(Ew$Ew[Ew_ind], na.rm=TRUE)
    }
    
    # inital concentrations
    cinit= c(M1 = df$M1[1], 
             M2 = df$M2[1], 
             M3 = df$M3[1], 
             M4 = df$M4[1], 
             M5 = df$M5[1])
    # time points for which conc is reported
    # include the points where data is available
    t=c(seq(min(df$time),max(df$time),1),df$time)
    t=sort(unique(t))
    # solve ODE for a given set of parameters
    out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
    
    # Filter data that contains time points where data is available
    outdf=data.frame(out)
    outdf=outdf[outdf$time %in% df$time,]
    # Evaluate predicted vs experimental residual
    preddf= outdf %>% 
      pivot_longer(-time, values_to = "conc", names_to = "species") %>% 
      arrange(time, species)
    
    expdf = df %>% 
      pivot_longer(-time, values_to = "conc", names_to = "species") %>% 
      arrange(time, species)
    
    ssqres=preddf$conc-expdf$conc
    SSQRES = c(SSQRES, ssqres)
  }
  
  # return predicted vs experimental residual
  return(SSQRES)
  
}

# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
parms=list(
  t1=0.15, 
  t2=0.15, 
  t3=0.15, 
  t4=0.15, 
  t5=0.15,
  m1=0.05, 
  m2=0.05, 
  m3=0.05, 
  m4=0.05, 
  m5=0.2)
# fitting
fitval=nls.lm(
  par=parms,
  lower=rep(0, length(parms)),
  upper=rep(1, length(parms)),
  fn=ssq)

# fitting
summary(fitval)
