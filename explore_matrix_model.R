library(tidyverse)
library(lubridate)
library(deSolve)

# mite life history
d1 = "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet=1)

lifestagenames = c("Ew", paste0('M', 1:5))

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
  filter(time>220) %>% #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mutate(Ew = ifelse(is.na(Ew), 0, Ew)) %>%
  identity()

# here I am trying to just use the previous stage to infer demographic parameters. I tried using ode with real eggs as the forcing parameter but its didnt do great - see "parameter_fitting.R" (perhaps because there are "impossible situation" where there are more of a lifestage than possible based on previous timepoints... measurement error?)

# derivatives for for reference
# dM1 = t1 * Ew - t2 * M1 - m1*M1
# dM2 = t2 * M1 - t3 * M2 - m2*M2 
# dM3 = t3 * M2 - t4 * M3 - m3*M3
# dM4 = t4 * M3 - t5 * M4 - m4*M4
# dM5 = t5 * M4 - m5*M5


dlag = d %>% 
  mutate(diff = time - lag(time, 1)) %>% 
  mutate(
    Ewlag = lag(Ew, 1),
    M1lag = lag(M1, 1),
    M2lag = lag(M2, 1),
    M3lag = lag(M3, 1), 
    M4lag = lag(M4, 1),
    M5lag = lag(M5, 1)
  ) %>% 
  filter(diff <= 14) 

# Ew 
mod0 = lm(Ew ~ Ewlag + M5lag - 1, data=dlag, weights = log(M5lag + 1))
summary(mod0)
car::Anova(mod0)
ggplot(dlag) +
  # partial residual plot
  geom_point(aes(M5lag, Ew - Ewlag*coef(mod0)['Ewlag'] )  ) + 
  geom_abline(slope=coef(mod0)['M5lag'], intercept = 0) 


# M1 
mod1 = lm(M1 ~ Ewlag + M1lag - 1, data=dlag, weights = log(M1+1))
summary(mod1)
car::Anova(mod1)
dlag %>% 
  filter(Ewlag != 0 )%>%
  ggplot() +
  # partial residual plot
  geom_point(aes(Ewlag, M1 - coef(mod1)['M1lag']*M1lag)) + 
  geom_line(aes(Ewlag, coef(mod1)['Ewlag']*Ewlag )) 


# M2 
mod2 = lm(M2 ~ Ewlag + M1lag + M2lag - 1, data=dlag, weights = log(M2+1))
summary(mod2)
car::Anova(mod2)
dlag %>% 
  filter(M1lag != 0 )%>%
  ggplot() +
  # partial residual plot
  geom_point(aes(M1lag, M2 - 
                   M2lag*coef(mod2)['M2lag'] - 
                   Ewlag*coef(mod2)['Ewlag'])) + 
  geom_line(aes(M1lag, coef(mod2)['M1lag'] * M1lag))


# M3 
mod3 = lm(M3 ~ M1lag + M2lag + M3lag - 1, data=dlag, weights = log(M3+1))
summary(mod3)
car::Anova(mod3)
dlag %>% 
  filter(M2lag != 0 )%>%
  ggplot() +
  # partial residual plot
  geom_point(aes(M2lag, M3 - 
                   M3lag*coef(mod3)['M3lag'] - 
                   M1lag*coef(mod3)['M1lag'])) + 
  geom_line(aes(M2lag, coef(mod3)['M2lag'] * M2lag)) 

# M4 
mod4 = lm(M4 ~ M2lag + M3lag + M4lag - 1, data=dlag, weights = log(M4+1))
summary(mod4)
car::Anova(mod4)
dlag %>% 
  filter(M3lag != 0 )%>%
  ggplot() +
  # partial residual plot
  geom_point(aes(M3lag, M4 - 
                   M4lag*coef(mod4)['M4lag'] - 
                   M2lag*coef(mod4)['M2lag'])) + 
  geom_line(aes(M3lag, coef(mod4)['M3lag'] * M3lag))

# M5 
mod5 = lm(M5 ~ M3lag + M4lag + M5lag - 1, data=dlag, weights = log(M5+1))
summary(mod5)
car::Anova(mod5)
dlag %>% 
  filter(M4lag != 0 )%>%
  ggplot() +
  # partial residual plot
  geom_point(aes(M4lag, M5 - 
                   M5lag*coef(mod5)['M5lag'] -
                   M3lag*coef(mod5)['M3lag'], color=time)) + 
  geom_line(aes(M4lag, coef(mod5)['M4lag'] * M4lag))  +
  scale_color_gradientn(colors=c('red', 'orange', 'yellow', 'green', 'blue'))



############### Build matrix model ##################
leslie = bind_rows(
  coef(mod0),
  coef(mod1),
  coef(mod2),
  coef(mod3),
  coef(mod4),
  coef(mod5),
) %>%
  as.matrix()
lesnames = colnames(leslie)[order(colnames(leslie))]
leslie = leslie[, lesnames]
leslie[is.na(leslie)] = 0
colnames(leslie) = c('Ew', paste0('M', 1:5))
rownames(leslie) = c('Ew', paste0('M', 1:5))
leslie

# simple simulation
tmax = 30
sim = matrix(NA, ncol=6, nrow=tmax)
sim[1, ] = c(0, 0, 0, 0, 0, 1000) 
for (i in 2:nrow(sim)) {
  sim[i, ] = leslie %*% sim[i-1, ]
}
sim %>% 
  as_tibble() %>%
  mutate(time = 1:tmax) %>% 
  pivot_longer(-time) %>% 
  ggplot() +
  geom_line(aes(time, value, color=name))


# simulation to predict data
id=4
di = d %>% 
  filter(group_id == id) 
dl = di %>%
  pivot_longer(-c(year:group_id)) 

time = seq(min(dl$time), max(dl$time), by=7)
sim = matrix(NA, ncol=6, nrow=length(time))
colnames(sim) = lifestagenames
sim[1, ] = as.numeric(di[1, lifestagenames])
forcing_var = 'Ew'
forcing  = dl %>% 
  filter(name==forcing_var) %>% 
  pull(value)
sim[,forcing_var] = forcing[1:length(time)] #ignores some data!!! 
col_update = !lifestagenames %in% forcing_var
for (i in 2:nrow(sim)) {
  sim[i, col_update] = (leslie %*% sim[i-1, ])[col_update]
}
siml = sim %>% 
  as_tibble() %>%
  mutate(time = time) %>% 
  pivot_longer(-time)

date0 = as.Date("2020-01-01") - 1
ggplot() +
  geom_line(data=siml, aes(date0+time, value, color=name)) +
  geom_point(data=dl, aes(date0+time, value, color=name)) +
  facet_wrap(~name)


