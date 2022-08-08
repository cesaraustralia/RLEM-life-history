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

