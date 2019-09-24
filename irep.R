# setwd('~/Documents/p_t_meta/')


library(tidyverse)
library(broom)


dat <- read_tsv('ALL_RESULTS_BULK_MAP.txt', skip = 1, col_types = c('iccnnnnnnninn'))%>%
  mutate(sample = sub('_mapped.sam','',sample), 
         genome = sub('.fa','',genome),
         day = sub('d([0-9][0-9])([a-z]+)([0-9][0-9][0-9])','\\1',sample), 
         treatment = sub('d([0-9][0-9])([a-z]+)([0-9][0-9][0-9])','\\2',sample), 
         bird = sub('d([0-9][0-9])([a-z]+)([0-9][0-9][0-9])','\\3',sample))  %>% 
  filter(!(is.na(iRep))) %>%
  filter(!grepl('ex', sample)) %>% 
  select(-X1, -`fragments/Mbp`)

# 3439 valid iRep values with 3x min cov... maybe should bump back up to 5
# 3342 actually without the 'extra' samples

dat5 <- dat %>% filter(coverage > 5)

# 2337 valid values with 5x cutoff
# 2278 
dat %>% ggplot(aes(x=log2(coverage))) + geom_histogram(bins = 100)

valid_irep_by_day <- dat %>% group_by(genome, day) %>% tally() %>% spread(key = day, value=n)
valid_irep_totbin <- dat %>% group_by(genome) %>% tally()


#
# this is not parsed correctly #
colnams <- c('bin', 'marker_lineage', 'num_genomes', 'num_markers', 'num_marker_sets', 'x0', 'x1', 'x2', 'x3','x4','x5+', 'Completeness', 'Contamination', 'Strain_heterogeneity')
checkm <- read_delim('CheckM_clean3.txt', delim = '\t', trim_ws = TRUE, skip = 2, col_names = colnams)


nrow(dat)

c(0,NA,4,6) >4


#### These ones we can probably do group comparisons at each day
look <- dat %>% 
  group_by(day, treatment, genome) %>% 
  tally() %>% spread(key = treatment, value=n) %>%
  ungroup() %>% group_by(day, genome) %>% 
  mutate(num_NA = sum(is.na(c(ctrl, sub, ther))), 
         num_w_4p = sum(c(ctrl, sub, ther) >3, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite(col = 'bin_day', genome, day)
  
group_compare <- look %>% filter(num_NA < 2) %>% filter(num_w_4p > 1)  

time_course <- look %>% filter(num_w_4p == 3)

# treatments by day
# only includes bins that have at least 2 treatments with 4+ observations
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% group_compare$bin_day) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day)

# days by treatment
# only includes bins with 4+ observations in all 3 treatments
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% time_course$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment)


dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% group_compare$bin_day) %>%
  group_by(genome, day, treatment) %>% 
  summarise(iRep=mean(iRep)) %>% 
  ggplot(aes(x=day, y=iRep, group=genome, color=treatment)) + geom_point() + geom_line() + 
  facet_wrap(~treatment)

 
dat %>%
  unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% group_compare$bin_day) %>%
  group_by(day, treatment, genome) %>%
  summarise(mean_irep=mean(iRep)) %>% 
  arrange(desc(mean_irep)) %>% 
  ggplot(aes(y=genome, x=mean_irep, color=treatment)) + geom_point() + facet_wrap(~day)


dat %>%
  unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% group_compare$bin_day) %>%
  group_by(day, treatment, genome) %>%
  summarise(mean_irep=mean(iRep)) %>% 
  arrange(desc(mean_irep)) %>% 
  filter(day =='07') %>% 
  ggplot(aes(y=genome, x=mean_irep, color=treatment)) + geom_point()#+ facet_wrap(~day)

#### need to start splitting these out by bin....


# install.packages("pillar")
# install.packages('dplyr')
# 
# nesty_dat <- dat %>%
#   unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
#   filter(bin_day %in% group_compare$bin_day) %>% 
#   filter(day == '07') %>% group_by(genome) %>%
#   nest()

nesty_dat <- dat %>%
  unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% group_compare$bin_day) %>% 
  group_by(genome, day) %>%
  nest()




test <- nesty_dat[25,][[3]][[1]]

summary(aov(data = test, formula = iRep ~ treatment))


tests_treat <- nesty_dat %>% 
  mutate(inter = map(data, ~ aov(data = ., formula = iRep ~ treatment)), 
         summ = map(inter, summary), 
         tid_sum = map(inter, tidy)) %>% 
  select(genome, tid_sum) %>% unnest(cols = c(tid_sum)) %>% 
  filter(term == 'treatment')


tests_treat %>% select(genome)




sessionInfo()



mtcars %>%
  split(.$cyl) %>%
  map(~ lm(mpg ~ wt, data = .)) %>%
  map(summary) %>%
  map_dbl("r.squared")


# install.packages('stringr')
# install.packages('tidyr')
# install.packages('purrr')
# install.packages('readr')
# install.packages('dplyr')
# install.packages('tibble')
# install.packages('tidyverse')
