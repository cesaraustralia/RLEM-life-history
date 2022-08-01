library(tidyverse)
library(lubridate)
library(deSolve)

# mite life history
d1 = "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet=1)

d = d1 %>% 
  mutate(DATE = as.Date(DATE)) %>%
  mutate(DATE = as.Date(format(DATE, "1990-%m-%d"))) %>%
  mutate(time = yday(DATE)) %>%
  select(c(-`ADULT MALE`, -`ADULT FEMALE`)) %>% 
  rename(
    year = YEAR,
    site = SITE, 
    block = BLOCK,
    sample = `SAMPLE No`,
    Ew = `EGGS ON PASTURE`,
    M1 = LARVA, 
    M2 = PROTO, 
    M3 = DEUTO, 
    M4 = TRITO,
    M5 = `TOTAL ADULTS`
  ) %>%  
  group_by(year, site, block) %>% 
  mutate(group_id = cur_group_id()) %>% 
  ungroup %>% 
  select(year, site, block, time, group_id, Ew, M1, M2, M3, M4, M5) %>%
  group_by(year, site, block, time, group_id) %>% 
  summarise_all(.funs = function(x) mean(x, na.rm=TRUE)) %>%
  ungroup %>%
  # filter(group_id == 1) %>% 
  mutate(Ew = ifelse(is.na(Ew), 0, Ew)) %>%
  identity()

# My plan is to keep substituting  in the previous terms so that the coefficients are more interpretable

# the assumption is that we can predict the current number of a stage class from a linear combination of the lagged stage classes and the lagged previous stage class. 
# for e.g. currentM2 = previousM2 - deathM2 - M2_to_M3 + M1_to_M2

# currentM2 = M2(t)
# previousM2 = M2(t-1)
# deathM2 = d*M2(t-1)
# M2_to_M3 = t3*M2(t-1)
# M1_to_M2 = t2*M1(t-1)

# Thus,
# M5 = M5lag - d*M5lag + t5*M4lag
# M5 = (1 - d)*M5lag + t5*M4lag

# M4 = (1- d)*M4lag + t4*M3lag + t5*M4lag  
# M3 = (1- d)*M3lag + t3*M2lag + t4*M3lag  
# M2 = (1- d)*M2lag + t2*M1lag + t3*M2lag  
# M1 = (1- d)*M1lag + t1*Ewlag + t2*M1lag  

# t4*M3lag = M4 - (1- d)*M4lag - t5*M4lag  
# t3*M2lag = M3 - (1- d)*M3lag - t4*M3lag  
# t2*M1lag = M2 - (1- d)*M2lag - t3*M2lag  
# t1*Ewlag = M1 - (1- d)*M1lag - t2*M1lag  


# M1 = (1- d)*M1lag + t1*Ewlag + M2 - (1- d)*M2lag - M3 + (1- d)*M3lag + M4 - (1- d)*M4lag - M5 + (1 - d)*M5lag   


# M2 
mod1.1 = lm(M1 ~ M1lag + Ewlag  + 
              offset(M2) + M2lag + 
              offset(M3) + M3lag + 
              offset(M4) + M4lag + 
              offset(M5) + M5lag - 1, data=dlag, weights = log(M2+1))
summary(mod1.1)
car::Anova(mod2)
ggplot(dlag) +
  # partial residual plot
  geom_point(aes(M1lag, M2 - M2lag*coef(mod2)['M2lag'])) + 
  geom_abline(slope=coef(mod2)['M1lag'], intercept = 0) +  
  scale_x_log10() +
  scale_y_log10()





d %>% 
  # ungroup %>%
  ggplot() + 
  geom_point(aes(M2, lag(M1, 1))) + 
  scale_x_log10() +
  scale_y_log10() + 
  geom_abline(slopr=1, intercept = 0)


 dl = d %>% 
  pivot_longer(c(-time, -year, -site, -block, -group_id)) 

# model


get_Ew = function(t) {
  Ew_ind = which(abs(d$time - t) == min(abs(d$time - t)))
  mean(d$Ew[Ew_ind], na.rm=TRUE)
}


rlem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
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


# eventfunc <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {
#       X=100
#     c(X, Y, Z)
#   })
# }


params <- c(
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

# ssq = function(p){
#   params=p
  init <- c(M1 = 0, M2 = 0, M3 = 0, M4 = 0, M5 = 0)
  tim <- seq(min(dl$time), max(dl$time), by = 1)
  sol <- ode(y = init, times = tim, func = rlem, 
             # events=list(func=eventfunc, time=10),
             parms = params 
             ) %>% 
    as.matrix() %>%
    as_tibble %>%
    mutate(time = as.numeric(time)) %>% 
    pivot_longer(-time, values_to = 'pred') %>% 
    mutate(pred = as.numeric(pred))
  
  # sum of squared deviations
  out = dl %>% 
    left_join(sol, by = c('time', 'name')) %>% 
    filter(name != "Ew") %>%
    mutate(resid = log(value+1) - log(pred+1)) %>% 
    drop_na(resid) %>%
    pull(resid) 
  return(out)

  ggplot() +
  geom_point(data=dl, aes(time, value, color=name)) +
  geom_line(data = sol, aes(time, pred, color=name)) + 
  facet_wrap(~name)

# }
# initial guess for parameters
# fitting
library(minpack.lm)
fitval=nls.lm(par=params,
              lower=rep(0, length(params)),
              upper=rep(1, length(params)),
              fn=ssq ,
              control = nls.lm.control(factor=100, maxiter=100, ptol=1e-9))
summary(fitval)

# 
# m1 =  nls(value~ssq(params), 
#       start=params)
# 
# 
# plot(sol[,'Y'], sol[,'X'], ty = 'l')
# plot(sol, which = c('Y','Z'))
# sol.lsoda <- lsoda(y = vars, times = tim, func = Lorenz, parms = params)
# plot(sol.lsoda, which = 'X')
