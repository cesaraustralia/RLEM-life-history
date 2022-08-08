library(tidyverse)
library(lubridate)
library(gam)

r <- read_csv("data/Ruesink.csv")

# Developmental time of stage j,
stages <- 1:2
Tj <- c(10, 8)
#  n must be made relatively large to obtain meaningful survival esti- mates. Increasing n has the effect of smoothing the f_S(t) curve
stepsize <- 10
df <- 9 # df for gam spline
tmin <- 10
tmax <- 130
mod <- list()
for (i in stages) {
  # fit simple gam
  mod[[i]] <- lm(formula(sprintf("N%d ~ bs(t ,df=%d)", i, df)), data = r)
}

# = Population density of stage j (var) at time t,
f_N <- function(i, t) {
  y <- pmax(0, predict(mod[[i]], newdata = data.frame(t = t)))
  y[t < tmin] <- 0
  y[t > tmax] <- 0
  return(y)
}
f_N(1, 10)
# Total number of individuals that have entered stage j up to time t, equation 1

f_C <- function(i, tt) {
  f_c <- function(i, t) {
    if (t >= Tj[i]) {
      f_N(i, t) + f_c(i, t - Tj[i])
    } else {
      f_N(i, t)
    }
  }
  sapply(tt, f_c, i = i)
}
f_C(1, seq(10, 40, by = 10))

# Total number of individuals that have left stage j up to time t, and
f_D <- function(i, t) {
  f_C(i, t - Tj[i])
}
# Rate individuals are entering stage j at time t
f_Q <- function(i, t) {
  n <- stepsize
  (f_C(i, t + n) - f_C(i, t - n)) / (2 * n)

}

#  Rate individuals are leaving stage j at time t
f_R <- function(i, t) {
  n <- stepsize
  (f_D(i, t + n) - f_D(i, t - n)) / (2 * n) 
}

# Survival between stages j and j+1
f_S <- function(i, t) {
  f_Q(i + 1, t) / f_R(i, t)
}

f_S(1, 20)


pred <- tibble(
  t = 10:130,
  N1 = f_N(1, t),
  N2 = f_N(2, t),
  S12 = f_S(1, t)* 100
  )  %>%
  pivot_longer(-t)

r %>%
  pivot_longer(-t) %>%
  ggplot(aes(t, value, color = name)) +
  geom_point() + 
  geom_line(data=pred) + 
  ylim(0, 100)

# recreate table 1 in Ruesink 1975
table1 <- tibble(
  t = r$t,
  N1 = f_N(1, t),
  C1 = f_C(1, t),
  D1 = f_D(1, t),
  R1 = f_R(1, t),
  S12 = f_S(1, t),
  Q2 = f_Q(2, t),
  C2 = f_C(2, t),
  D2 = f_D(2, t),
  N2 = f_N(2, t)
)

table1




