
###### Purpose of this script is to identify genomic bins that:
#               1) Have different growth rates between treatments
#               2) Have different growth rates over time
#



library(tidyverse)
library(broom)


dat <- read_tsv('ALL_RESULTS_BULK_MAP.txt', skip = 1, col_types = c('iccnnnnnnninn'))%>%
  mutate(sample = sub('_mapped.sam','',sample), 
         genome = sub('.fa','',genome),
         day = sub('d([0-9][0-9])([a-z]+)([0-9][0-9][0-9])','\\1',sample), 
         treatment = sub('d([0-9][0-9])([a-z]+)([0-9][0-9][0-9])','\\2',sample), 
         bird = sub('d([0-9][0-9])([a-z]+)([0-9][0-9][0-9])','\\3',sample), 
         day_num = as.numeric(day))  %>% 
  filter(!(is.na(iRep))) %>%
  filter(!grepl('ex', sample)) %>% 
  select(-X1, -`fragments/Mbp`)


# 3439 valid iRep values with 3x min cov... maybe should bump back up to 5
# 3342 actually without the 'extra' samples

dat5 <- dat %>% filter(coverage > 5)
dat %>% ggplot(aes(x=log2(coverage))) + geom_histogram(bins = 100)


# 2337 valid values with 5x cutoff
# 2278 


#### iRep is significantly negatively correlated with both coverage and relative abundance


cor.test(x = dat$iRep, log(dat$`relative abundance`))
dat %>% ggplot(aes(x=log(`relative abundance`), y=iRep)) + geom_point() + geom_smooth(method = 'lm')


cor.test(x = dat$iRep, log(dat$coverage))
dat %>% ggplot(aes(x=log(coverage), y=iRep)) + geom_point() + geom_smooth(method = 'lm')
log(5)


###### relabund and coverage v. strongly correlated #####
cor.test(x = log(dat$`relative abundance`), log(dat$coverage))
dat %>% ggplot(aes(x=log(coverage), y=log(`relative abundance`))) + geom_point() + geom_smooth(method = 'lm')







valid_irep_by_day <- dat %>% group_by(genome, day) %>% tally() %>% spread(key = day, value=n)
valid_irep_totbin <- dat %>% group_by(genome) %>% tally() %>% arrange(desc(n))


#  CheckM bin info  #
colnams <- c('bin', 'marker_lineage', 'num_genomes', 'num_markers', 'num_marker_sets', 'x0', 'x1', 'x2', 'x3','x4','x5+', 'Completeness', 'Contamination', 'Strain_heterogeneity')
checkm <- read_delim('CheckM_clean3.txt', delim = '\t', trim_ws = TRUE, skip = 2, col_names = colnams)


# calculating which bins have enough datapoints available for meaningful stats

look <- dat %>% 
  group_by(day, treatment, genome) %>% 
  tally() %>% spread(key = treatment, value=n) %>%
  ungroup() %>% group_by(day, genome) %>% 
  mutate(num_NA = sum(is.na(c(ctrl, sub, ther))), 
         num_w_4p = sum(c(ctrl, sub, ther) >3, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite(col = 'bin_day', genome, day, remove = FALSE)
  
#### These have at least 2 groups with 4+ observations
two_groups <- look %>%  filter(num_w_4p > 1)  

# these have all 3 groups represented by at least 4 data points each
all_groups <- look %>% filter(num_w_4p == 3)

# no bin has 4 or more observations in every group at every timepoint
look %>% select(genome, day, num_w_4p) %>%
  spread(key = day, value = num_w_4p) %>% 
  mutate(nine = `07`+ `35`+ `78`) %>% arrange(desc(nine)) %>% 
  mutate(d7g = ifelse(`07` == 3, TRUE, FALSE), 
         d35g = ifelse(`35` == 3, TRUE, FALSE),
         d78g = ifelse(`78` == 3, TRUE, FALSE))


## one timepoint ##

# 16 bins have all treatment groups represented with 4+ observations at D7
D7 <- all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  filter(!is.na(`07`))


# 22 bins have all treatment groups represented with 4+ observations at D35
D35 <- all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  filter(!is.na(`35`))

# 16 bins have all treatment groups represented with 4+ observations at D78
D78 <- all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  filter(!is.na(`78`))


### 2 timepoints ###

# 10 bins have all 3 groups with 4+ observations at both d7 and d35
all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  mutate(d7vd35 = `07` + `35`) %>% 
  filter(!is.na(d7vd35))

# 2 bins have all 3 groups with 4+ observations at both d35 and d78
all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  mutate(d35vd78 = `35` + `78`) %>% 
  filter(!is.na(d35vd78))


### 4+ in 2 groups by day ###


# these ones can compare d7 vs d35 in ctrl

#
ctrl735 <- look %>% 
#  filter(num_w_4p >1) %>%
#  filter(ctrl >3) %>%
  select(-bin_day, -num_NA, -num_w_4p) %>% 
  gather(key = 'treat', value = 'count', -day, -genome) %>%
  unite(col='treat_day', treat, day) %>% spread(key=treat_day, value=count) %>% 
  select(genome, starts_with('ctrl')) %>% 
  filter(ctrl_07 > 3 & ctrl_35 > 3)

sub735 <- look %>% 
  #  filter(num_w_4p >1) %>%
  #  filter(ctrl >3) %>%
  select(-bin_day, -num_NA, -num_w_4p) %>% 
  gather(key = 'treat', value = 'count', -day, -genome) %>%
  unite(col='treat_day', treat, day) %>% spread(key=treat_day, value=count) %>% 
  select(genome, starts_with('sub')) %>% 
  filter(sub_07 > 3 & sub_35 > 3)

ther735 <- look %>% 
  #  filter(num_w_4p >1) %>%
  #  filter(ctrl >3) %>%
  select(-bin_day, -num_NA, -num_w_4p) %>% 
  gather(key = 'treat', value = 'count', -day, -genome) %>%
  unite(col='treat_day', treat, day) %>% spread(key=treat_day, value=count) %>% 
  select(genome, starts_with('ther')) %>% 
  filter(ther_07 > 3 & ther_35 > 3)


# ctrl735$genome

# these bins can be used to compare growth rates in all 3 treatments at days 7 and 35
lmbins <- intersect(intersect(ctrl735$genome, sub735$genome), ther735$genome)





d735_lms <- dat %>% filter(genome %in% lmbins & day %in%c('07','35')) %>% group_by(genome) %>% nest()

test <- d735_lms[4,2][[1]][[1]]
test %>% ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day)


summary(lm(data = test, formula = iRep ~ treatment*day_num))
########

nesty_dat <- dat %>% filter(genome %in% D7$genome & day %in%c('07') & treatment != 'sub') %>%
  group_by(genome) %>% nest() 


mayb <- nesty_dat %>% 
  mutate(inter = map(data, ~ pairwise.t.test(x = .$iRep, g= .$treatment, p.adjust.method = 'none', pool.sd = TRUE)), 
         #summ = map(inter, summary), 
         tid_sum = map(inter, tidy)) %>% 
  select(genome, tid_sum) %>% unnest(cols = c(tid_sum)) #%>% 
  #filter(term == 'treatment') #%>% filter(adj.p.value < 0.05)


mayb$p.adj <- p.adjust(mayb$p.value, method = 'fdr')


res <- mayb %>% filter(p.adj < 0.1)

checkm$genome <- sub('_new','',checkm$bin)


checkm[checkm$genome %in% res$genome,]



dat %>% filter(genome %in% res$genome[1]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[2]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[3]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[4]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[5]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[6]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[7]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[8]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[9]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()

dat %>% filter(genome %in% res$genome[10]) %>% filter(day == '07') %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + geom_point(aes(color=`relative abundance`), size=3) + scale_color_viridis_c()



test <- nesty_dat[2,2][[1]][[1]]
summary(pairwise.t.test(x=test$iRep, g=test$treatment))
summary(aov(data=test, formula=iRep ~ treatment))


#############

dat %>%
  mutate(day=as.numeric(day)) %>% 
  filter(genome == 'bin.716' & treatment =='ctrl' & day %in% c(7, 35)) %>%
  ggplot(aes(x=day, y=iRep)) +
  geom_point() +
  geom_smooth(method = 'lm')


test <- dat %>%
  mutate(day=as.numeric(day)) %>% 
  filter(genome == 'bin.716' & treatment =='ctrl' & day %in% c(7, 35))



summary(lm(data=test, formula = iRep ~day))


############ plots ###########

# treatments by day
# only includes bins that have at least 2 treatments with 4+ observations
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day)

# days by treatment
# only includes bins with 4+ observations in all 3 treatments
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment)


dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>%
  group_by(genome, day, treatment) %>% 
  summarise(iRep=mean(iRep)) %>% 
  ggplot(aes(x=day, y=iRep, group=genome, color=treatment)) + geom_point() + geom_line() + 
  facet_wrap(~treatment)

 
dat %>%
  unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>%
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

dat %>%
  unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% group_compare$bin_day) %>%
  group_by(day, treatment, genome) %>%
  summarise(mean_irep=mean(iRep)) %>% 
  arrange(desc(mean_irep)) %>% 
  filter(day =='78') %>% 
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

tests_treat$fdr <- p.adjust(tests_treat$p.value, method = 'fdr')
tests_treat$BH <- p.adjust(tests_treat$p.value, method = 'BH')
tests_treat$BY <- p.adjust(tests_treat$p.value, method = 'BY')
tests_treat$holm <- p.adjust(tests_treat$p.value, method = 'holm')
tests_treat$hoch <- p.adjust(tests_treat$p.value, method = 'hochberg')
tests_treat$homm <- p.adjust(tests_treat$p.value, method = 'hommel')
tests_treat$bonferroni <- p.adjust(tests_treat$p.value, method = 'bonferroni')


tests_treat %>% filter(BH < 0.1)

###
p.adjust.methods

dat


dat7 <- dat %>% filter(day == '07')

summary(lm(data = dat7, formula = iRep ~ genome*treatment))




# 
# sessionInfo()
# 
# 
# 

# mtcars %>%
#   split(.$cyl) %>%
#   map(~ lm(mpg ~ wt, data = .)) %>%
#   map(summary) %>%
#   map_dbl("r.squared")


# install.packages('stringr')
# install.packages('tidyr')
# install.packages('purrr')
# install.packages('readr')
# install.packages('dplyr')
# install.packages('tibble')
# install.packages('tidyverse')
