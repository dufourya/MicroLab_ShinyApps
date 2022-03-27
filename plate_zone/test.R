library(rjags)
library(ggplot2)
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())
library(tidybayes)

DF <-
  data.frame(
    Zone = c(0,10,20,10,0,15,25,11,2,8,22,12),
    Disk = c("A", "B", "C", "D","A", "B", "C", "D","A", "B", "C", "D"),
    Plate = c("1","1","1","1","2","2","2","2","3","3","3","3"),
    Student = c("JJ","JJ","JJ","JJ","JJ","JJ","JJ","JJ","CC","CC","CC","CC"),
    stringsAsFactors = FALSE
  )


DF$Disk = as.factor(DF$Disk)
DF$Plate = as.factor(DF$Plate)
DF$Student = as.factor(DF$Student)

data = compose_data(DF)

print(data)

params_inits = list(
  "sigma" = 10,
  "sigma_d" = 10,
  "sigma_p" = 10,
  "sigma_s" = 10
)

model <-
  jags.model(
    "plate_zone.bug",
    data = data,
    inits = params_inits,
    n.chains = 4
  )

update(model, n.iter = 10000)

samples <-
  coda.samples(
    model,
    variable.names = c(
      "d","p","s"
    ),
    n.iter = 2500,
    thin = 10
  )

samples = recover_types(samples,DF)

samples %>%
  spread_draws(d[Disk]) %>%
  ggplot(aes(x = Disk, y = d)) +
  geom_violin(color = NA, fill = "gray65") +
  stat_pointinterval(.width = c(.95, .66))

samples %>%
  spread_draws(p[Plate]) %>%
  ggplot(aes(x = Plate, y = p)) +
  geom_violin(color = NA, fill = "gray65") +
  stat_pointinterval(.width = c(.95, .66))

samples %>%
  spread_draws(s[Student]) %>%
  ggplot(aes(x = Student, y = s)) +
  geom_violin(color = NA, fill = "gray65") +
  stat_pointinterval(.width = c(.95, .66))


samples %>%
  spread_draws(d[Disk]) %>%
  compare_levels(d, by = Disk) %>%
  ggplot(aes(y = Disk, x = d)) + geom_eyeh() + geom_vline(xintercept = 0)

samples %>%
  spread_draws(p[Plate]) %>%
  compare_levels(p, by = Plate) %>%
  ggplot(aes(y = Plate, x = p)) + geom_eyeh() + geom_vline(xintercept = 0)

samples %>%
  spread_draws(s[Student]) %>%
  compare_levels(s, by = Student) %>%
  ggplot(aes(y = Student, x = s)) + geom_eyeh() + geom_vline(xintercept = 0)
