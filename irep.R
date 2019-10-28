
###### Purpose of this script is to identify genomic bins that:
#               1) Have different growth rates between treatments
#               2) Have different growth rates over time
#



library(tidyverse)
library(broom)


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


# dat %>% ggplot(aes(y=iRep, x=0)) +
#   geom_violin() +
#   geom_jitter(alpha=.05) +
#   ylim(1,4) + xlim(-1, 1) +
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank(), 
#         axis.title.x = element_blank())


dat %>% ggplot(aes(x=iRep)) +
  geom_histogram(bins = 50) +
  theme_bw()



dat %>% ggplot(aes(x=(coverage))) + geom_histogram(bins = 100)  # some bins have extremely high coverage
dat %>% ggplot(aes(x=log2(coverage))) + geom_histogram(bins = 100) 


# 2278 valid iRep estimates 


#### iRep is significantly negatively correlated with both coverage and relative abundance


dat$log_relabund <- log(dat$`relative abundance`)
dat %>% ggplot(aes(x=log_relabund, y=iRep)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

cor.test(x = dat$iRep, log(dat$`relative abundance`))



# valid_irep_by_day <- dat %>% group_by(genome, day) %>% tally() %>% spread(key = day, value=n)
valid_irep_totbin <- dat %>% group_by(genome) %>% tally() %>% arrange(desc(n))

valid_irep_totbin %>%
  ggplot(aes(x=n)) +
  geom_histogram(bins = 50) +
  xlab('number of valid iRep estimates')+
  ggtitle('number of valid iRep estimates per bin') + 
  theme_bw()

#### NOW ONE WITH BINS THAT OCCUR 4+ TIMES IN EITHER 2 TIMEPOINTS PER TREAT or 4+ TIMS IN 2 TREATS PER TIMEPOINT



#####

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



### 4+ in 2 groups by day ###


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

# ctrl735$genome

# these bins can be used to compare growth rates in all 3 treatments at days 7 and 35
#4 bins
# probably a little too complicated to try and assess the effect of treatment as well as timepoint at 
# the same time, especialyl becuase the data are sparse....
lmbins <- intersect(intersect(ctrl735$genome, sub735$genome), ther735$genome)
lmbins


#### this is looking
test <- dat %>% filter(genome %in% lmbins & day %in%c('07','35'))

summary(lm(data = test, formula = iRep ~ treatment*day))

summary(lm(data = test, formula = iRep ~ genome + day + treatment))

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


dat %>% filter(genome %in% lmbins & day %in%c('07','35')) %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment, shape=day)) + geom_boxplot() +geom_jitter()



dat %>% filter(genome == 'bin.493' & day %in%c('07','35')) %>% 
  ggplot(aes(x=day, y=iRep, fill=treatment), shape=21) + geom_boxplot() +geom_jitter(position = position_dodge2(width = .75))




# bin.493: no terms seem to influence growth rate.  no evidence that treatment group or day 

###########
### CONTROL ONLY DIFFERENCES BETWEEN GROWTH RATES AT D7 AND D35

dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=day)) + geom_boxplot() 


dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35'))%>% 
  ggplot(aes(x=day, y=iRep, fill=day)) + geom_violin() 


test <- dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35'))

summary(lm(data=test, formula = iRep ~ day+genome))

ctrl735_lms <- dat %>% filter(genome %in% ctrl735$genome & day %in% c('07', '35')) %>% group_by(genome) %>% nest()

ctrl735_lms %>% mutate(lms=map(data, ~ lm(data=. , formula = iRep ~ day)), 
                       tid_sum = map(lms, tidy)) %>% select(genome, tid_sum) %>% 
  unnest(cols = c('tid_sum')) %>% filter(term == 'day35' & p.value < 0.05)

############ DAY 7 DIFFS BTWEEN GRUOPS ##########

dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  ggplot(aes(x=treatment, y=iRep, fill=treatment)) + geom_violin() + geom_jitter(shape=21, width = .2)+
  ggtitle('iRep growth rate estimates at D7', '') + theme_bw() +
  stat_summary(fun.y = "mean", colour = "black", size = 2, geom = "point")

dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) + 
  ggtitle('iRep growth rate estimates at D7') + theme_bw()


dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  group_by(genome, treatment) %>% 
  summarise(miRep=mean(iRep), 
            se=sd(iRep)/sqrt(n())) %>% 
  ggplot(aes(x=genome, y=miRep, fill=treatment)) +
  geom_col(position = 'dodge', color='black') + 
  geom_errorbar(aes(ymin=miRep-se, ymax=miRep + se), position = 'dodge')


D7$genome

D7_treat_comp <- dat %>% filter(genome %in% D7$genome & day %in% c('07'))


sigs <- D7_treat_comp %>% group_by(genome) %>% nest() %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ treatment)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) %>% 
  select(-adj.p.value) %>% filter(tuk_pval < 0.05)

sigs

library(lme4)
library(lmerTest)





dat %>% filter(genome %in% D7$genome & day %in% c('07'))%>% 
  filter(genome %in% sigs$genome) %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) + 
  ggtitle('iRep growth rate estimates at D7') + theme_bw()


summary(lmer(data = D7_treat_comp, formula = iRep ~ treatment + (1+treatment|genome)))
summary(lmer(data = D7_treat_comp, formula = iRep ~ treatment + (1|genome)))



summary(lm(data=D7_treat_comp, formula = iRep ~ treatment + genome))
summary(lm(data=D7_treat_comp, formula = iRep ~ treatment * genome))
summary(lm(data=D7_treat_comp, formula = iRep ~ genome * treatment))

checkm %>% filter(bin %in% sigs$genome)


###### DAY 35 DIFF BTWEEN GROUPS #######

D35$genome

D35_treat_comp <- dat %>% filter(genome %in% D35$genome & day %in% c('35'))


summary(lm(data=D35_treat_comp, formula = iRep ~ treatment + genome))
summary(lm(data=D35_treat_comp, formula = iRep ~ treatment * genome))


D35_treat_comp %>% ggplot(aes(x=treatment, y=iRep, fill=treatment)) +
  geom_violin() + geom_jitter(shape=21, width = .2) + 
  ggtitle('iRep growth rate estimates at D35') + theme_bw()


dat %>% filter(genome %in% D35$genome & day %in% c('35'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) +
  ggtitle('iRep growth rate estimates at D35, by genome') + theme_bw()




sigs <- D35_treat_comp %>% group_by(genome) %>% nest() %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ treatment)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) %>% 
  select(-adj.p.value) %>% filter(tuk_pval < 0.05)

sigs



###### Day 78 Diff BTWEEN GROUPS ######


dat %>% filter(genome %in% D78$genome & day %in% c('78'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() + geom_jitter(shape=21, position=position_dodge2(width = .75))

D78$genome

D78_treat_comp <- dat %>% filter(genome %in% D78$genome & day %in% c('78'))


summary(lm(data=D78_treat_comp, formula = iRep ~ treatment + genome))
summary(lm(data=D78_treat_comp, formula = iRep ~ treatment * genome))


D78_treat_comp %>% ggplot(aes(x=treatment, y=iRep, fill=treatment)) +
  geom_violin() + geom_jitter(shape=21, width = .2) + 
  ggtitle('iRep growth rate estimates at D78') + theme_bw()


D78_treat_comp %>% #filter(genome %in% D35$genome & day %in% c('35'))%>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() +
  geom_jitter(shape=21, position=position_dodge2(width = .75)) +
  ggtitle('iRep growth rate estimates at D78, by genome') + theme_bw()


sigs <- D78_treat_comp %>% group_by(genome) %>% nest() %>% 
  mutate(ANOVA=map(data, ~ aov(data=., iRep ~ treatment)), 
         tuk = map(ANOVA, TukeyHSD), 
         tid_tuk = map(tuk, tidy)) %>% 
  select(genome, tid_tuk) %>% 
  unnest(cols = tid_tuk) %>%
  mutate(tuk_pval = adj.p.value, 
         fdr.pval = p.adjust(tuk_pval, method = 'fdr')) %>% 
  select(-adj.p.value) %>% filter(tuk_pval < 0.05)


sigs











########  Excluding subther?  #######

nesty_dat <- dat %>% filter(genome %in% D7$genome & day %in%c('07') & treatment != 'sub') %>%
  group_by(genome) %>% nest() 

nesty_dat <- dat %>% filter(genome %in% D35$genome & day %in%c('35') & treatment != 'sub') %>%
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

### D7 boxplots some sig diff btw groups

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


##### D35 boxplots

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

### all data, unfiltered ###

# histogram
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  ggplot(aes(x=iRep)) + 
  geom_histogram(bins=50)+
  theme_bw()

# boxplots
dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day)

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment)


dat %>% nrow()
# 2278 observations
dat %>% select(genome) %>% unlist() %>% unique() %>% length()
# 244 genomes



### 1 with 4 plus ###
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
  ggtitle('NEED TITLE',
          'considering only those genomes that have \nat least 1 treatment with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% one_group$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment) +
  ggtitle('NEED TITLE', 
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
  ggtitle('NEED TITLE',
          'considering only those genomes that have \nat least 2 treatments with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment) +
  ggtitle('NEED TITLE', 
          'considering only those genomes that have \nat least 2 treatments with 4+ observations at any time')


#### DELETE BELOW HERE TO
#HIST

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=treatment, y=iRep)) + geom_boxplot() + facet_wrap(~day)

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% two_groups$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment)


### HERE


### 3 with 4 plus ###
# only includes bins with 4+ observations in all 3 treatments

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
  ggtitle('NEED TITLE',
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')

dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>% 
  ggplot(aes(x=day, y=iRep)) + geom_boxplot() + facet_wrap(~treatment) +
  ggtitle('NEED TITLE', 
          'considering only those genomes that have \nat least 3 treatments with 4+ observations at any time')



###########


dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  filter(bin_day %in% all_groups$bin_day) %>%group_by(genome) %>% 
  arrange(desc(iRep)) %>% ungroup() %>% 
  mutate(genome2=factor(genome, levels = unique(genome))) %>% 
  ggplot(aes(x=genome, y=iRep, fill=treatment)) + geom_boxplot() + geom_jitter(alpha=.2, shape=21)


dat %>% unite(col = 'bin_day', genome, day, remove = FALSE) %>% 
  # filter(bin_day %in% all_groups$bin_day) %>%
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

