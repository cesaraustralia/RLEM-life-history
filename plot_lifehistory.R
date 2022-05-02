library(tidyverse)
library(lubridate)

# mite life history
d1 = "data/mite population_ ngn_keys_1990-1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet=1)
names(d1)
lifestagenames = c("LARVA", "PROTO", "DEUTO", "TRITO", "ADULT FEMALE", 
                   "ADULT MALE", "TOTAL ADULTS", "EGGS ON PASTURE")
d1 %>% 
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


# mite eggs
d2 = "data/lifehistory/egg development_ ngn_keys_1990_1992 20_Dec.xlsx" %>%
  readxl::read_xlsx(sheet=1)
d2$SITE %>% unique

d2 %>% 
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


# over summering eggs
d3 = "data/summer sample_ngn_keys_1990-1992 17-Jan.xlsx" %>%
  readxl::read_xlsx(sheet=1) %>% 
  mutate(`No. EGGS` = as.numeric(`No. EGGS`) ) %>% 
  mutate(DATE= as.Date(DATE))

dsum = d3 %>% 
  group_by(SITE, SAMPLE, DATE) %>% 
  summarise(nbodies = sum(`BODIES/CORE`,na.rm=T),
            nbodiesdissected = n(), 
            total_eggs = sum(`No. EGGS`, na.rm=T),
            mean_eggs_per_mite = mean(`No. EGGS`, na.rm=T), 
            se_eggs = sd(`No. EGGS`)/sqrt(nbodiesdissected),
            .groups="drop") %>% 
  mutate(nbodiesdissected = if_else(nbodies == 0, as.integer(0), nbodiesdissected))

dsum %>% 
  pivot_longer(cols = c(nbodies, nbodiesdissected, total_eggs, mean_eggs_per_mite, se_eggs)) %>% 
  filter(name %in% c("nbodies", "mean_eggs_per_mite")) %>%
  filter(value < 400) %>% # !!! filter outlier for plotting
  mutate(name = recode(name, 
                       nbodies = "number of mites", 
                       mean_eggs_per_mite = "mean eggs per mite")) %>% 
  mutate(Year = factor(year(DATE))) %>%
  ggplot(aes(Year, value, group=DATE)) +
  geom_violin() +
  # scale_y_log10() +
  facet_grid(name~SITE, switch = "y", scales="free_y")

ggsave("plots/lifehistoryegg_summer.png", width = 10, height=7)


