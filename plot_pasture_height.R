library(tidyverse)
library(lubridate)

# effect of pasture height
col_types = 
  c("date","numeric","text", rep("numeric", 34), 
   "text", "text", "numeric", "numeric", "text")  

d = "data/Short Tall mites and pasture 17-Jan.xls" %>%
  readxl::read_xls(sheet=1, 
                   col_types = col_types) %>% 
  mutate(DATE= as.Date(DATE)) %>%
  mutate(CLOV_CM_SHORT = as.numeric(CLOV_CM_SHORT)) %>%
  mutate(CLOV_CM_TALL = as.numeric(CLOV_CM_TALL)) %>%
  pivot_longer(cols = -c(DATE, WEEK, SAMPLE, YEAR,	
                         State,	SITE,	Latitude,	Longitude,	COLLECTOR)) %>% 
  drop_na(value) %>% 
  mutate(name = gsub("TALL", "_TALL", name)) %>%
  mutate(name = gsub("SHORT", "_SHORT", name)) %>%
  separate(col = name, into=c("var", "pasture"), sep = "__")

names(d)

dsum = d %>% 
  group_by(SITE, DATE, WEEK, pasture, var) %>% 
  summarise(
    mean_val = mean(value), 
    se_val = sd(value)/sqrt(n()),
    .groups="drop") 

unique(d$var)
dsum %>% 
  ggplot(aes(DATE, mean_val, col=pasture)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_val-se_val, ymax=mean_val+se_val, width=1)) + 
  geom_line() +
  facet_wrap(~var, scales = "free_y") + 
  theme_bw() + 
  xlab("Date") + 
  ylab("Mean value")

ggsave("plots/pasture_height_effect.png", width = 4, height=3, scale = 2)  
