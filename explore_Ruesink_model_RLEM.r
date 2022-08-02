# mite life history
source('explore_Ruesink_model.R')
library(tidyverse)


d1 <- "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet = 1)

lifestagenames <- c("Ew", paste0("M", 1:5))

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
    N1 = `EGGS ON PASTURE`,
    N2 = LARVA,
    N3 = PROTO,
    N4 = DEUTO,
    N5 = TRITO,
    N6 = `TOTAL ADULTS`
  ) %>%
  group_by(year, site) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  select(year, site, block, time, group_id, N1, N2, N3, N4, N5, N6) %>%
  group_by(year, site, time, group_id) %>%
  summarise_all(.funs = function(x) mean(x, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(group_id == 2) %>%
  # filter(time>220) %>% #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mutate(N1 = ifelse(is.na(N1), 0, N1)) %>%
  identity()

r <- d %>% 
select(t = time, N1, N2, N3, N4, N5, N6)

# Developmental time of stage j,
stages <- 1:6
Tj <- c(5, 5, 5, 5, 5, 5)
#  n must be made relatively large to obtain meaningful survival esti- mates. Increasing n has the effect of smoothing the f_S(t) curve
stepsize <- 7
df <- 9 # df for gam spline
tmin <- min(r$t)
tmax <- max(r$t)
mod <- list()
for (i in stages) {
  # fit simple gam
  mod[[i]] <- lm(formula(sprintf("N%d ~ bs(t ,df=%d)", i, df)), data = r)
}


pred <- tibble(
  t = tmin:tmax,
  N1 = f_N(1, t),
  N2 = f_N(2, t),
  N3 = f_N(3, t),
  N4 = f_N(4, t),
  N5 = f_N(5, t),
  N6 = f_N(6, t),
#   S12 = f_S(1, t)* 100
  )  %>%
  pivot_longer(-t)

r %>%
  pivot_longer(-t) %>%
  ggplot(aes(t, value, color = name, linetype=name)) + 
  scale_color_viridis_d(end=0.8, direction=-1) + 
#   geom_point() +
  geom_line(data=pred)  
  facet_wrap(~name)
