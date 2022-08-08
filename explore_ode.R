library(tidyverse)
library(lubridate)
library(deSolve)

#### OPTIONS ##########
# comment out options if using explore_ode_get_pars.R
# id <- 1 # site number 

# mite life history
narr <- read_csv("data/narrogin_19210101to20211231.csv")
keys <- read_csv("data/keysbrook_19210101to20211231.csv")
clim <- bind_rows(
  mutate(narr, SITE = "KEYSBROOK (DELBORELLO)"),
  mutate(keys, SITE = "NARROGIN (CURNOW)")
) %>%
  mutate(DATE = as.Date(`YYYY-MM-DD`, format = "%Y-%m-%d")) %>%
  mutate(YEAR = year(DATE)) %>%
  filter(YEAR > 1989) %>%
  mutate(TEMP = (min_temp + max_temp) / 2) %>%
  group_by(SITE, YEAR) %>%
  mutate(dd = cumsum(TEMP)) %>%
  ungroup() %>%
  select(SITE, DATE, dd)

d1 <- "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet = 1)

d <- d1 %>%
  mutate(DATE = as.Date(DATE)) %>%
  left_join(clim, by = c("DATE", "SITE")) %>%
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
  group_by(year, site) %>% # group block
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  select(year, site, time = dd, group_id, Ew, M1, M2, M3, M4, M5) %>%
  group_by(year, site, time, group_id) %>%
  summarise_all(.funs = function(x) mean(x, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Ew = ifelse(is.na(Ew), 0, Ew)) %>%
  identity()


rlem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # dEd = 0
    # dEw = eggrate * M5 / 2 - t1 * Ew
    Ew <- get_Ew(t)
    # Ew = 100
    dM1 <- t1 * Ew - t2 * M1 - m1 * M1
    dM2 <- t2 * M1 - t3 * M2 - m2 * M2
    dM3 <- t3 * M2 - t4 * M3 - m3 * M3
    dM4 <- t4 * M3 - t5 * M4 - m4 * M4
    dM5 <- t5 * M4 - m5 * M5
    dEw <-  e * M5 - t1 * Ew - m0 * Ew 
    list(c(dM1, dM2, dM3, dM4, dM5, dEw))
  })
}

dd <- d %>%
  filter(group_id == id)

dl <- dd %>%
  select(time, Ew, M1, M2, M3, M4, M5) %>%
  pivot_longer(-time)


dEw <- dl %>%
  filter(name == "Ew") %>%
  select(time, value) %>%
  drop_na()
get_Ew <- function(t) approxfun(dEw$time, dEw$value, yright = 0.1, yleft = 0.1)(t)

dM5 <- dl %>%
  filter(name == "M5") %>%
  select(time, value) %>%
  drop_na()
get_M5 <- function(t) approxfun(dM5$time, dM5$value, yright = 0.1, yleft = 0.1)(t)


# when you use degree days the interpretation of these parameters changes. y = exp(0.05*t)
params <- c(
  t1 = 0.03,
  t2 = 0.005,
  t3 = 0.010,
  t4 = 0.015,
  t5 = 0.020,
  m1 = 0.01,
  m2 = 0.01,
  m3 = 0.01,
  m4 = 0.01,
  m5 = 0.01,
  m0 = 0.01,
  e  = 10
)

max_temp = 30
interpret_params= function(params, daily_temp){
   daily_temp * params
}
interpret_params(params, 15)

run_sim <- function(params) {
  init <- c(
    M1 = dd$M1[1], M2 = dd$M2[1], M3 = dd$M3[1],
    M4 = dd$M4[1], M5 = dd$M5[1], Ew = dd$Ew[1]
  )
  tim <- c(seq(min(dl$time), max(dl$time), by = 10), dd$time)
  tim <- tim[order(tim)]
  ode(
    y = init, times = tim, func = rlem,
    # events=list(func=eventfunc, time=10),
    parms = params
  ) %>%
    as.matrix() %>%
    as_tibble() %>%
    mutate(time = as.numeric(time)) %>%
    pivot_longer(-time, values_to = "pred") %>%
    mutate(pred = as.numeric(pred))
}
system.time(run_sim(params) )
########### PLOT MODEL ###########################

interpret_params(params, 15)

sol <- run_sim(params) %>% 
  # filter(name!="Ew") %>% 
  mutate(value = pred)

dl %>%
  # left_join(sol, by=c('time','name')) %>%
  ggplot(aes(time, value, color = name)) +
  geom_point() +
  geom_line(data=sol) +
  # scale_y_log10() +
  facet_wrap(~name, scales = "free")

sol$pred %>% log() %>% sum
log(dl$value + 1) %>% sum

############ OPTIMISE PARAMETERS ##############
# make objective function
obj <- function(params) {
  cat(".")
  print(params)
  sol <- run_sim(params) %>% 
    mutate(pred = pmax(0, pred))
 

  resid = dl %>%
    # filter(name == "M5") %>% # only consider adults
    # weight by instar
    mutate(weights = as.numeric(str_extract(name, "\\d+"))) %>%
    mutate(weights = ifelse(is.na(weights), 0, weights)) %>%
    mutate(weights = (time - min(time))  / max(time - min(time)) * weights) %>%
    left_join(sol, by = c("time", "name")) %>%
    mutate(resid = (log(value + 1) - log(pred + 1))^2) %>%
    # mutate(resid = (value - pred)^2) %>%
    mutate(resid = weights * resid) %>%
    pull(resid)

    sum(resid)
}




res <- optim(
  params,
  obj,
  method = "L-BFGS-B",
  # Bound pars at 0 to 30 (100% dev or transition at 30C)
  lower = c(rep(0.0001, 5), rep(-1/max_temp, 5),    0,       1),
  upper = c(rep(1/max_temp, 5), rep(1/max_temp, 5), 0.0001, 10),
  control = list(maxit = 2)
)
params <- res$par

######## COPY PASTE PLOT CODE FROM ABOVE ########
sol <- run_sim(params) %>% 
  # filter(name!="Ew") %>% 
  mutate(value = pred)

p1 = dl %>%
  # left_join(sol, by=c('time','name')) %>%
  ggplot(aes(time, value, color = name)) +
  geom_point() +
  geom_line(data=filter(sol,name!="Ew")) +
  # scale_y_log10() +
  facet_wrap(~name)
p1
ggsave(sprintf('plots/sim_fit/site_id_%s.png', id), p1)
