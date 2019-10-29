library(tidyverse)
library(broom)
library(lme4)
library(lmerTest)

#  CheckM bin info  #
colnams <- c('bin', 'marker_lineage', 'num_genomes', 'num_markers', 'num_marker_sets', 'x0', 'x1', 'x2', 'x3','x4','x5+', 'Completeness', 'Contamination', 'Strain_heterogeneity')

checkm <- read_delim('CheckM_clean3.txt', delim = '\t', trim_ws = TRUE, skip = 2, col_names = colnams) %>%
  mutate(bin=sub('_new', '', bin)) 



# read in the data and extract some metadata from sample names
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


dat <- dat %>% filter(coverage > 5)  # only keep estimates from samples where bins had > 5x coverage (as reccommended)

# Histogram of all iRep estimates
dat %>% ggplot(aes(x=iRep)) +
  geom_histogram(bins = 50) +
  theme_bw() + ggtitle('Histogram of all valid iRep estimates')


# some bins have extremely high coverage / relative abundance
dat %>% ggplot(aes(x=(coverage))) + geom_histogram(bins = 100) + 
  ggtitle('Histogram of bin coverages') +
  theme_bw()

dat %>%
  ggplot(aes(x=log2(coverage))) +
  geom_histogram(bins = 100) +
  ggtitle('Histogram of log(coverage)') + 
  theme_bw()


dat %>%
  ggplot(aes(x=`relative abundance`)) +
  geom_histogram(bins = 100) +
  ggtitle('Histogram of relative_abundance') + 
  theme_bw()

# dat$log_relabund <- log(dat$`relative abundance`)

dat %>%
  ggplot(aes(x=(log2(`relative abundance`)))) +
  geom_histogram(bins = 100) +
  ggtitle('Histogram of log(relative_abundance)') + 
  theme_bw()

#### iRep is significantly negatively correlated with both coverage and relative abundance

cor.test(x = dat$iRep, log(dat$`relative abundance`))
dat %>% ggplot(aes(x=log(`relative abundance`), y=iRep)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

# what about testing this within each bin?


valid_irep_totbin <- dat %>% 
  group_by(genome) %>%
  tally() %>%
  arrange(desc(n))

valid_irep_totbin %>%
  ggplot(aes(x=n)) +
  geom_histogram(bins = 50) +
  xlab('number of valid iRep estimates')+
  ggtitle('number of valid iRep estimates per bin') + 
  theme_bw()


dat %>% nrow()
# 2278 valid iRep estimates

dat %>% select(genome) %>% unlist() %>% unique() %>% length()
# from 244 genomes


# How do these estimates look by treatment? by day?
# boxplots
dat %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day) + 
  theme_bw()

dat %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment)+ 
  theme_bw()

### neat, maybe a treatment effect? 
### maybe a time effect? 

### But, we dont know if we are measuring the same genomes growth rates in each of these categories...

#####  main point here is that the data are sparse #######

# calculating which bins have enough datapoints available for meaningful stats

number_of_obs <- dat %>% 
  group_by(day, treatment, genome) %>% 
  tally() %>% spread(key = treatment, value=n) %>%
  ungroup() %>% group_by(day, genome) %>% 
  mutate(num_NA = sum(is.na(c(ctrl, sub, ther))), 
         num_w_4p = sum(c(ctrl, sub, ther) >3, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite(col = 'bin_day', genome, day, remove = FALSE)
  
number_of_obs


# data grouped by bin and day to tally numobs per group


# These have at least 1 treatment group with 4+ observations at the given timepoint

one_group <- number_of_obs %>%  filter(num_w_4p >= 1)
one_group


#### These have at least 2 groups with 4+ observations
two_groups <- number_of_obs %>%  filter(num_w_4p > 1)  
two_groups

# these have all 3 groups represented by at least 4 data points each
all_groups <- number_of_obs %>% filter(num_w_4p == 3)
all_groups

# no bin has 4 or more observations in every group at every timepoint
number_of_obs %>% select(genome, day, num_w_4p) %>%
  spread(key = day, value = num_w_4p) %>% 
  mutate(all_times = `07`+ `35`+ `78`) %>% arrange(desc(all_times)) %>% 
  mutate(d7g = ifelse(`07` == 3, TRUE, FALSE), 
         d35g = ifelse(`35` == 3, TRUE, FALSE),
         d78g = ifelse(`78` == 3, TRUE, FALSE))


# kindof makes sense, both ABX treatment and time affected the community composition
#
# Will probably have to stick to comparing iRep estimates in two main ways:
#  1) Within timepoint between treatments
#  2) within a treatment between timepoints


## one timepoint ##

# 8 bins have all treatment groups represented with 4+ observations at D7
D7 <- all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  filter(!is.na(`07`))
D7

# 11 bins have all treatment groups represented with 4+ observations at D35
D35 <- all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  filter(!is.na(`35`))
D35

# 8 bins have all treatment groups represented with 4+ observations at D78
D78 <- all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  filter(!is.na(`78`))
D78

### 2 timepoints ###

# 4 bins have all 3 groups with 4+ observations at both d7 and d35
all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  mutate(d7vd35 = `07` + `35`) %>% 
  filter(!is.na(d7vd35))

# 1 bin has all 3 groups with 4+ observations at both d35 and d78
all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  mutate(d35vd78 = `35` + `78`) %>% 
  filter(!is.na(d35vd78))

# no bin has all 3 groups with 4+ observations at d7 and d78
all_groups %>%
  select(genome, day, num_w_4p) %>% 
  spread(key=day, value = num_w_4p) %>% 
  mutate(d07vd78 = `07` + `78`) %>% 
  filter(!is.na(d07vd78))



##### These sets of bins can be used to investigate time effect #####
##### within each treatment #####

# these ones can compare d7 vs d35 in ctrl

# 15 bins
ctrl735 <- number_of_obs %>% 
#  filter(num_w_4p >1) %>%
#  filter(ctrl >3) %>%
  select(-bin_day, -num_NA, -num_w_4p) %>% 
  gather(key = 'treat', value = 'count', -day, -genome) %>%
  unite(col='treat_day', treat, day) %>% spread(key=treat_day, value=count) %>% 
  select(genome, starts_with('ctrl')) %>% 
  filter(ctrl_07 > 3 & ctrl_35 > 3)
ctrl735

# 8 bins
sub735 <- number_of_obs %>% 
  #  filter(num_w_4p >1) %>%
  #  filter(ctrl >3) %>%
  select(-bin_day, -num_NA, -num_w_4p) %>% 
  gather(key = 'treat', value = 'count', -day, -genome) %>%
  unite(col='treat_day', treat, day) %>% spread(key=treat_day, value=count) %>% 
  select(genome, starts_with('sub')) %>% 
  filter(sub_07 > 3 & sub_35 > 3)
sub735

# 11 bins
ther735 <- number_of_obs %>% 
  #  filter(num_w_4p >1) %>%
  #  filter(ctrl >3) %>%
  select(-bin_day, -num_NA, -num_w_4p) %>% 
  gather(key = 'treat', value = 'count', -day, -genome) %>%
  unite(col='treat_day', treat, day) %>% spread(key=treat_day, value=count) %>% 
  select(genome, starts_with('ther')) %>% 
  filter(ther_07 > 3 & ther_35 > 3)
ther735



# start to pair down the data to only include genomes that have enough data for comparisons?

##### only includes estimates from bins that have at least 1 treatment with 4+ observations in any one timepoint  ##


dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>% nrow()
# 1572 observations

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>%
  select(genome) %>% unlist() %>% unique() %>% length()
#102 genomes


# histogram
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>% 
  ggplot(aes(x=iRep)) + 
  geom_histogram(bins=50)+
  theme_bw() + 
  ggtitle('Histogram of iRep estimates',
          'considering only those genomes that have \nat least 1 treatment with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>%
  group_by(genome) %>%
  tally() %>% 
  ggplot(aes(x=n))+geom_histogram(bins=30) + 
  ggtitle('Valid iRep estimates per genome', 
          'considering only those genomes that have \nat least 1 treatment with 4+ observations at any time')


# boxplots
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day) + 
  ggtitle('iRep estimates by treatment',
          'considering only those genomes that have \nat least 1 treatment with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment) +
  ggtitle('iRep estimates by time', 
          'considering only those genomes that have \nat least 1 treatment with 4+ observations at any time')



### 2 with 4 plus ###
# only includes bins that have at least 2 treatments with 4+ observations

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% nrow()
# 1062 observations

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>%
  select(genome) %>% unlist() %>% unique() %>% length()
#56 genomes


# histogram
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=iRep)) + 
  geom_histogram(bins=50)+
  theme_bw() + 
  ggtitle('Histogram of iRep estimates',
          'considering only those genomes that have \nat least 2 treatments with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>%
  group_by(genome) %>%
  tally() %>% 
  ggplot(aes(x=n))+geom_histogram(bins=30) + 
  ggtitle('Valid iRep estimates per genome', 
          'considering only those genomes that have \nat least 2 treatments with 4+ observations at any time')


# boxplots
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day) + 
  ggtitle('iRep estimates by treatment',
          'considering only those genomes that have \nat least 2 treatments with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment) +
  ggtitle('iRep estimates by time', 
          'considering only those genomes that have \nat least 2 treatments with 4+ observations at any time')



### 3 with 4 plus ###
# only includes bins with 4+ observations in all 3 treatments at any one timepoint

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>% nrow()
# 507 observations

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>%
  select(genome) %>% unlist() %>% unique() %>% length()
#22 genomes


# histogram
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>% 
  ggplot(aes(x=iRep)) + 
  geom_histogram(bins=50)+
  theme_bw() + 
  ggtitle('Histogram of iRep estimates',
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>%
  group_by(genome) %>%
  tally() %>% 
  ggplot(aes(x=n))+geom_histogram(bins=30) + 
  ggtitle('Valid iRep estimates per genome', 
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')


# boxplots
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day) + 
  ggtitle('iRep estimates by treatment',
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment) +
  ggtitle('iRep estimates by time', 
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')

### Well, the data aren't idea, but we'll try to subset it and get something meaningful out of it.

# specifically I am interested if there is statistical evidence for these two hypothesis:
# 1) growth rates are supressed in the abx treatments relative to the controls
# 2) growth rates are lower in the d35 communities relative to the D7 communities

####### TREATMENT EFFECTS ###########
# only trying to asses treatment effect within each day. 

############ DAY 7 DIFFS BTWEEN GRUOPS ##########
# only includes genomes with 4+ observations in all 3 treatment groups at D7
# 8 genomes, 149 observations

dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  ggplot(aes(x=treatment, y=iRep, fill=treatment)) + geom_violin() + geom_jitter(shape=21, width = .2)+
  ggtitle('iRep growth rate estimates at D7', '') + theme_bw() +
  stat_summary(fun.y = "mean", colour = "black", size = 2, geom = "point")

dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) + 
  ggtitle('iRep growth rate estimates at D7') + theme_bw()


D7_treat_comp <- dat %>% filter(genome %in% D7$genome & day %in% c('07'))


D7_tukeys <- D7_treat_comp %>% group_by(genome) %>% nest() %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ treatment)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr'))

sigs <- D7_tukeys %>%
  select(-adj.p.value) %>% filter(tuk_pval < 0.05)

sigs
sigs_fdr <- sigs %>% filter(fdr.pval < 0.1)


# only those bins with a sig difference detected in tukey's test
dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  filter(genome %in% sigs$genome) %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) + 
  ggtitle('iRep growth rate estimates at D7', 'tukey pval < 0.05 -- not adjusted for multiple tukeys') + theme_bw()



dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  filter(genome %in% sigs_fdr$genome) %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) + 
  ggtitle('iRep growth rate estimates at D7', 'fdr pval < 0.10') + theme_bw()





D7_tukeys %>% #filter(comparison != 'ther-sub') %>% 
  ggplot(aes(x=genome, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2) + 
  coord_flip() + geom_hline(yintercept = 0, color='red') + 
  ggtitle('95% confidence intervals for the difference in growth rate between treatments',
          'negative estimates indicate higher growth rate in ctrl, D7 only')+ facet_wrap(~comparison)

# mixed model? not super confident about this.
summary(lmer(data = D7_treat_comp, formula = iRep ~ treatment + (1|genome)))
checkm %>% filter(bin %in% sigs$genome) %>% select(bin, marker_lineage, Completeness, Contamination)




###### DAY 35 DIFF BTWEEN GROUPS #######
# only includes genomes with 4+ observations in all 3 treatment groups at D35
# 11 genomes, 211 observations


D35_treat_comp <- dat %>% filter(genome %in% D35$genome & day %in% c('35'))


D35_treat_comp %>% ggplot(aes(x=treatment, y=iRep, fill=treatment)) +
  geom_violin() + geom_jitter(shape=21, width = .2) + 
  ggtitle('iRep growth rate estimates at D35') + theme_bw()


dat %>% filter(genome %in% D35$genome & day %in% c('35'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) +
  ggtitle('iRep growth rate estimates at D35, by genome') + theme_bw()




D35_tukeys <- D35_treat_comp %>% group_by(genome) %>% nest() %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ treatment)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) 

sigs <- D35_tukeys %>% 
  select(-adj.p.value) %>% filter(tuk_pval < 0.05)

# no sigs by anova/tukey


D35_tukeys %>% #filter(comparison != 'ther-sub') %>% 
  ggplot(aes(x=genome, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2) + 
  coord_flip() + geom_hline(yintercept = 0, color='red') + 
  ggtitle('95% confidence intervals for the difference in growth rate between treatments',
          'negative estimates indicate higher growth rate in ctrl, D35 only')+ facet_wrap(~comparison)

# mixed model? not super confident about this.
summary(lmer(data = D35_treat_comp, formula = iRep ~ treatment + (1|genome)))

# no treatment effect detectable at D35




###### Day 78 Diff BTWEEN GROUPS ######
# only includes genomes with 4+ observations in all 3 treatment groups at D35
# 8 genomes, 147 observations

dat %>% filter(genome %in% D78$genome & day %in% c('78'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() + geom_jitter(shape=21, position=position_dodge2(width = .75))


D78_treat_comp <- dat %>% filter(genome %in% D78$genome & day %in% c('78'))



D78_treat_comp %>% ggplot(aes(x=treatment, y=iRep, fill=treatment)) +
  geom_violin() + geom_jitter(shape=21, width = .2) + 
  ggtitle('iRep growth rate estimates at D78') + theme_bw()


D78_treat_comp %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) +
  ggtitle('iRep growth rate estimates at D78, by genome') + theme_bw()


D78_tukeys <- D78_treat_comp %>% group_by(genome) %>% nest() %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ treatment)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr'))


sigs <- D78_tukeys %>% 
  select(-adj.p.value) %>% filter(tuk_pval < 0.05)


D78_tukeys %>% filter(comparison != 'ther-sub') %>% 
  ggplot(aes(x=genome, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2) + 
  coord_flip() + geom_hline(yintercept = 0, color='red') + 
  ggtitle('95% confidence intervals for the difference in growth rate between treatments',
          'negative estimates indicate higher growth rate in ctrl, D35 only')+ facet_wrap(~comparison)



# mixed model? not super confident about this.
summary(lmer(data = D78_treat_comp, formula = iRep ~ treatment + (1|genome)))
# small treatment effect detected, both sub and ther, but at this timepoint they are receiving equal doses.

########### END TREATMENT EFFECT #########



####### BEGIN TIME EFFECT ##########



### CONTROL ONLY DIFFERENCES BETWEEN GROWTH RATES AT D7 AND D35
# only includes genomes with 4+ observations at both D7 and D35 in the control treatment only
# 15 genomes 201 observations

dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35')& treatment == 'ctrl')%>% 
  ggplot(aes(x=day, y=iRep, fill=day)) +
  geom_violin() +
  ggtitle('iRep estimates at D7 and D35, control only', '15 genomes, 201 observations')

# looks to be a trend to higher growth rates in d7 communities

dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35') & treatment == 'ctrl')%>% 
  ggplot(aes(x=genome, y=iRep, fill=day)) + geom_boxplot() +
  ggtitle('iRep estimates at D7 and D35, control only', '15 genomes, 201 observations')

ctrl735_dat <- dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35')& treatment == 'ctrl')

ctrl735_tests <- dat %>% 
  filter(genome %in% ctrl735$genome & day %in% c('07', '35')& treatment == 'ctrl') %>%
  group_by(genome) %>%
  nest()


# fit a simple linear model on each genome
# ctrl735_lms <- ctrl735_tests %>% mutate(lms=map(data, ~ lm(data=. , formula = iRep ~ day)), 
#                        tid_sum = map(lms, tidy)) %>% select(genome, tid_sum) %>% 
#   unnest(cols = c('tid_sum'))#%>% filter(term == 'day35' & p.value < 0.05)
# 
# ctrl735_lms %>%
#   filter(term == 'day35' & p.value < 0.05)

# Of the 15 genomes, there is evidence that 5 of them have lower growth rates at d35 relative to D7


# ANOVA for each genome

ctrl735_tuk <- ctrl735_tests %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ day)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) %>% 
  select(-adj.p.value) #%>% filter(tuk_pval < 0.05)

ctrl735_tuk %>% 
  filter(tuk_pval < 0.05)


ctrl735_tuk %>%
  ggplot(aes(x=genome, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2) + 
  coord_flip() + geom_hline(yintercept = 0, color='red') + 
  ggtitle('95% confidence intervals for the difference in growth rate between D7 and D35',
          'negative estimates indicate higher growth rate at D7, Control treatment only')

summary(lmer(data = ctrl735_dat, formula = iRep ~ day + (1|genome)))

# Time effect detectable in the control group.  Lower growth-rate at D35 compared to D7 


### SUBTHER ONLY DIFFERENCES BETWEEN GROWTH RATES AT D7 AND D35
# only includes genomes with 4+ observations at both D7 and D35 in the control treatment only
# 8 genomes 94 observations
dat %>% filter(genome %in% sub735$genome & day %in% c('07', '35')& treatment == 'sub')%>% 
  ggplot(aes(x=genome, y=iRep, fill=day)) + geom_boxplot() +
  ggtitle('iRep estimates at D7 and D35, sub only', '8 genomes, 94 observations')


dat %>% filter(genome %in% sub735$genome & day %in% c('07', '35')& treatment == 'sub')%>% 
  ggplot(aes(x=day, y=iRep, fill=day)) + geom_violin() +
  ggtitle('iRep estimates at D7 and D35, sub only', '8 genomes, 94 observations')

# dat$treatment
sub735_dat <- dat %>% filter(genome %in% sub735$genome & day %in% c('07', '35')& treatment == 'sub')


sub735_tests <- dat %>% 
  filter(genome %in% sub735$genome & day %in% c('07', '35')& treatment == 'sub') %>%
  group_by(genome) %>%
  nest()


# ANOVA for each genome

sub735_tuk <- sub735_tests %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ day)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) %>% 
  select(-adj.p.value) #%>% filter(tuk_pval < 0.05)

sub735_tuk %>% 
  filter(tuk_pval < 0.05)


sub735_tuk %>%
  ggplot(aes(x=genome, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2) + 
  coord_flip() + geom_hline(yintercept = 0, color='red') + 
  ggtitle('95% confidence intervals for the difference in growth rate between D7 and D35',
          'negative estimates indicate higher growth rate at D7, Control treatment only')

summary(lmer(data = sub735_dat, formula = iRep ~ day + (1|genome)))

# Time effect detectable in the sub group.  Lower growth-rate at D35 compared to D7 


### THER ONLY DIFFERENCES BETWEEN GROWTH RATES AT D7 AND D35
# only includes genomes with 4+ observations at both D7 and D35 in the control treatment only
# 11 genomes 131 observations

dat %>% filter(genome %in% ther735$genome & day %in% c('07', '35')& treatment == 'ther')%>% 
  ggplot(aes(x=genome, y=iRep, fill=day)) + geom_boxplot() +
  ggtitle('iRep estimates at D7 and D35, ther only', '11 genomes, 131 observations')



dat %>% filter(genome %in% ther735$genome & day %in% c('07', '35')& treatment == 'ther')%>% 
  ggplot(aes(x=day, y=iRep, fill=day)) + geom_violin() +
  ggtitle('iRep estimates at D7 and D35, ther only', '11 genomes, 131 observations')


# dat$treatment
ther735_dat <- dat %>% filter(genome %in% ther735$genome & day %in% c('07', '35')& treatment == 'ther')


ther735_tests <- dat %>% 
  filter(genome %in% ther735$genome & day %in% c('07', '35')& treatment == 'ther') %>%
  group_by(genome) %>%
  nest()



# ANOVA for each genome

ther735_tuk <- ther735_tests %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ day)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) %>% 
  select(-adj.p.value) #%>% filter(tuk_pval < 0.05)

ther735_tuk %>% 
  filter(tuk_pval < 0.05)

# Of the 11 genomes, there is evidence that 0 of them have lower growth rates at d35 relative to D7



ther735_tuk %>%
  ggplot(aes(x=genome, y=estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2) + 
  coord_flip() + geom_hline(yintercept = 0, color='red') + 
  ggtitle('95% confidence intervals for the difference in growth rate between D7 and D35',
          'negative estimates indicate higher growth rate at D7, Therapeutic treatment only')

summary(lmer(data = ther735_dat, formula = iRep ~ day + (1|genome)))


# No time effect in ther








############ DUMP BELOW HERE #########
# these are all half baked (or unbaked) ideas,  dont look at them.


ctrl_735_lmod <- lmer(data = ctrl735_dat, formula = iRep ~ day + (1|genome))
summary(ctrl_735_lmod)
confint(ctrl_735_lmod)
lme4::ranef(ctrl_735_lmod)


sub_735_lmod <- lmer(data = sub735_dat, formula = iRep ~ day + (1|genome))
summary(sub_735_lmod)
confint(sub_735_lmod)
lme4::ranef(sub_735_lmod)


ther_735_lmod <- lmer(data = ther735_dat, formula = iRep ~ day + (1|genome))
summary(sub_735_lmod)
confint(sub_735_lmod)
lme4::ranef(sub_735_lmod)


######### From above ########

######## bin in either of these three groups #####
# get an idea of a d07 to d35 time effect?
# these are bins that have 4+ observations at both D7 and D35 within any treatment

all735 <- unique(c(ctrl735$genome, sub735$genome, ther735$genome))

dat_all735 <- dat %>% filter(genome %in% all735 & day %in%c('07','35'))

dat_all735 %>% ggplot(aes(x=day, y=iRep, fill=day)) +
  geom_boxplot() +
  stat_summary(geom = 'point', fun.y = 'mean', color='blue') + 
  facet_wrap(~treatment) + 
  theme_bw()

dat_all735 %>%
  ggplot(aes(x=genome, y=iRep, color=day)) +
  geom_boxplot() + coord_flip()+
  facet_wrap(~treatment, scales = 'free')

day_difs <- dat_all735 %>% group_by(genome, day, treatment) %>% 
  summarise(miRep=mean(iRep))

dat_all735 %>% group_by(genome, day) %>% 
  summarise(miRep=mean(iRep)) %>%
  ggplot(aes(x=genome, y=miRep, color=day)) + geom_point()+
  geom_segment(aes(yend=miRep, xend=genome))+
  coord_flip()

####### show direction of all d7 vs d35 changes

summary(lmer(data=dat_all735, formula = iRep ~ day+treatment + (1|genome)))

# summary(lmer(data = D7_treat_comp, formula = iRep ~ treatment + (1|genome)))





##### Probably going to limit the detection of a time effect to between D7 and D35

#### bin in all of these three groups ####
# ctrl735$genome

# these bins can be used to compare growth rates in all 3 treatments at days 7 and 35
#4 bins
# probably a little too complicated to try and assess the effect of treatment as well as timepoint at 
# the same time, especialyl becuase the data are sparse....
lmbins <- intersect(intersect(ctrl735$genome, sub735$genome), ther735$genome)
lmbins


#### this is looking at only those bins that have 4+ observations in all three treatment groups at days 7 and 35
# 4 genomes, 161 observations

lmbins_dat <- dat %>% filter(genome %in% lmbins & day %in%c('07','35'))


lmbins_dat %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() + facet_wrap(~day) +
  ggtitle('iRep estimates by time', 
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')


# summary(lm(data = lmbins_dat, formula = iRep ~ genome + treatment*day))
# summary(lm(data = lmbins_dat, formula = iRep ~ genome + day + treatment))

d735_lms <- dat %>% filter(genome %in% lmbins & day %in%c('07','35')) %>% group_by(genome) %>% nest()

# d735_lms[4,1]
# test <- d735_lms[3,2][[1]][[1]]
# test %>% ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day)



###########

hmmm <- d735_lms %>%ungroup() %>% 
  mutate(lms=map(data, ~ lm(data=. , formula = iRep ~ treatment*day)), 
         tid_sum = map(lms, tidy)) %>% select(genome, tid_sum) %>% 
  unnest(cols = c('tid_sum'))


hmmm %>% filter(p.value <= 0.05)


# dat %>% filter(genome %in% lmbins & day %in%c('07','35')) %>% 
#   ggplot(aes(x=genome, y=iRep, fill=treatment, shape=day)) + geom_boxplot() +geom_jitter()



# dat %>% filter(genome == 'bin.493' & day %in%c('07','35')) %>% 
#   ggplot(aes(x=day, y=iRep, fill=treatment), shape=21) + geom_boxplot() +geom_jitter(position = position_dodge2(width = .75))
# 
# 


# bin.493: no terms seem to influence growth rate.  no evidence that treatment group or day 





