library(tidyverse)
library(lubridate)
library(deSolve)

# to do 
# add physiological time

# mite life history
d1 <- "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet = 1)

d <- d1 %>%
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
  group_by(year, site) %>% # group block
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  select(year, site, time, group_id, Ew, M1, M2, M3, M4, M5) %>%
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
    list(c(dM1, dM2, dM3, dM4, dM5))
  })
}


id = 5

dd = d %>% 
  filter(group_id == id)

dl <- dd %>%
  select(time, Ew, M1, M2, M3, M4, M5) %>%
  pivot_longer(-time)


get_Ew <- function(t) {
  xy <- dl %>%
    filter(name == "Ew") %>%
    select(time, value) %>%
    drop_na()
  approxfun(xy$time, xy$value)(t)
}

params <- c(
  t1 = 0.05,
  t2 = 0.05,
  t3 = 0.10,
  t4 = 0.15,
  t5 = 0.20,
  m1 = 0.01,
  m2 = 0.01,
  m3 = 0.01,
  m4 = 0.01,
  m5 = 0.15
)


init <- c(M1 = 0, M2 = 0, M3 = 0, M4 = 0, M5 = 0)
tim <- seq(min(dl$time), max(dl$time), by = 7)
sol <- ode(
  y = init, times = tim, func = rlem,
  # events=list(func=eventfunc, time=10),
  parms = params
) %>%
  as.matrix() %>%
  as_tibble() %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(-time, values_to = "pred") %>%
  mutate(pred = as.numeric(pred))

ggplot() +
  geom_point(data = dl, aes(time, value, color = name)) +
  geom_line(data = sol, aes(time, pred, color = name)) +
  facet_wrap(~name)

