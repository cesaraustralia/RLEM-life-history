library(tidyverse)
library(lubridate)

# mite eggs
d1 = "data/lifehistory/egg development_ ngn_keys_1990_1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet=1)
d1$SITE %>% unique

d1 %>% 
  mutate(DATE = as.Date(DATE)) %>%
  mutate(diapause_prop = AESTIV./(WINTER + AESTIV.)) %>% 
  mutate(diapause = diapause_prop > 0.5) %>%
  group_by(SITE, YEAR, DATE) %>% 
  summarise(mean_diapause = mean(diapause) ,groups="drop") %>% 
  ggplot(aes(DATE, mean_diapause)) + 
  geom_point() + 
  geom_line() +
  scale_x_date(date_labels = "%b %Y") +
  facet_wrap( ~ SITE, scales = "free_x") 
ggsave("plots/lifehistoryegg.png", width = 10)


# mite life history
d2 = "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet=1)
names(d2)
lifestagenames = c("LARVA", "PROTO", "DEUTO", "TRITO", "ADULT FEMALE", 
                   "ADULT MALE", "TOTAL ADULTS", "EGGS ON PASTURE")
d2 %>% 
  mutate(DATE = as.Date(DATE)) %>%
  mutate(DATE2 = as.Date(format(DATE, "1990-%m-%d"))) %>%
  pivot_longer(cols = lifestagenames, names_to = "stage", values_to = "number") %>% 
  mutate(stage = factor(stage, levels = lifestagenames)) %>%
  filter(stage != "TOTAL ADULTS") %>%
  group_by(SITE, YEAR, DATE2, stage) %>% 
  summarise(mean_number = mean(number, na.rm=T)) %>% 
  ggplot(aes(DATE2, mean_number, fill=stage)) +
  geom_area() + 
  scale_fill_viridis_d() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  facet_grid(YEAR ~ SITE, scales = "free") +
  xlab("") +
  ylab("Abundance")
ggsave("plots/lifehistorymite.png", width = 10, height=7)

d = d2 %>% 
  drop_na(`EGGS ON PASTURE`) 

ggplot(d) +
  geom_point(aes(`ADULT FEMALE`, `EGGS ON PASTURE`)) + 
  scale_x_log10() + 
  scale_y_log10()

lm(data=d)
