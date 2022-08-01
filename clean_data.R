# clean data 
library(tidyverse)

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
  identity()
d
write_csv(d, 'data/clean_data.csv')