# get parameter estimates for Ridsdill-Smith 1997 data
library(tidyverse)

# loop over 5 site ids (Narr 2 years, Keys 3 years)

pars = tibble()
for (id in 1:5) {
    source('explore_ode.R')
    pars = pars %>% 
        bind_rows(
            as_tibble(t(params)) 
        ) 
       
}
temp_C = 15

write_csv(pars, "./plots/sim_fit/pars.csv")

psum = pars %>% 
    select(t1:m5) %>%
    mutate(no = 1:n()) %>%
    pivot_longer(-no) %>%
    select(-no) %>%
    # convert C.d^-1 to d^-1 at temp_C
    mutate(value = temp_C * value) %>% 
    group_by(name) %>%
    summarise_all(
        list(mu=mean, sd=sd, n=length)) %>% 
    mutate(se = sd/sqrt(n)) %>% 
    mutate(var = str_extract(name, '.')) %>% 
    mutate(name = str_extract(name, '(\\d+)'))


ggplot(psum, aes(mu, name, xmin=mu-sd, xmax=mu+sd)) + 
    geom_point() + 
    geom_errorbarh(height=0.1) + 
    facet_wrap(~var) + 
    ylab('RLEM instar') + 
    xlab('mean instart mortality and development rate (1/d)')

ggsave('plots/sim_fit/pars_mu_se.png', width=10, height=5)
