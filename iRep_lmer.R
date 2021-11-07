
library(tidyverse)
library(broom)
library(lme4)
library(lmerTest)
library(emmeans)




dat_QC <- read_tsv('iRep_estimates_QCd.tsv') |> 
  mutate(day=factor(day, levels = c('07', '35', '78')), 
         treatment=factor(treatment, levels = c('ctrl', 'sub', 'ther')))

dat_QC$day
dat_QC$treatment

mod_interact <- lmer(data=dat_QC, formula = iRep ~ treatment * day + (1|genome))
mod_simple <- lmer(data=dat_QC, formula = iRep ~ treatment + day + (1|genome))


summary(mod_simple)
summary(mod_interact)

# some of the interaction coefs are significant, ther group prob has different 
# behavior over time 

anova(mod_simple, mod_interact)
# interaction model fits the data better

plot(mod_simple)
plot(mod_interact)

resids <- mod_interact %>% resid()


LOOK <- dat_QC |> mutate(resids)

LOOK |> ggplot(aes(x=resids, y=iRep)) + geom_point()

treatment_effects <- emmeans(mod_interact, ~ treatment | day) %>%
  contrast(method='pairwise', adjust='fdr')

time_effects <- emmeans(mod_interact, ~ day | treatment) %>%
  contrast(method='pairwise', adjust='fdr') 


all_effects <- rbind(treatment_effects, time_effects) %>%tidy(conf.int=T) |> 
  mutate(contrast=factor(contrast))



library(cowplot)
# treatment effects at different timepoints
all_effects |>
  filter(treatment != '.') |>
  ggplot(aes(x=fct_rev(contrast), color=treatment, y=estimate, ymin=conf.low, ymax=conf.high)) + 
  geom_linerange(position=position_dodge(width = .51), size=1.25, ) + 
  geom_point(aes(fill=treatment),position=position_dodge(width = .51), shape=21, size=4, color='black')+
  coord_flip() + 
  geom_hline(yintercept = 0)+ 
  theme_cowplot() + 
  theme(panel.grid.major = element_line(color='grey')) + 
  xlab('contrast between days')+ 
  ylab('difference in iRep growth rates') + 
  ggtitle('Effects of time within each treatment') + 
  ylim(-.1, .3)

# time effects within different treatments
all_effects |>
  filter(day != '.') |>
  ggplot(aes(x=fct_rev(contrast),color=day, y=estimate, ymin=conf.low, ymax=conf.high)) + 
  geom_linerange(position=position_dodge(width = .51), size=1.25, ) + 
  geom_point(aes(fill=day),position=position_dodge(width = .51), shape=21, size=4, color='black')+
  coord_flip() + 
  geom_hline(yintercept = 0)+ 
  theme_cowplot() + 
  theme(panel.grid.major = element_line(color='grey'))+
  xlab('contrast between treatments')+ 
  ylab('difference in iRep growth rates') + 
  ggtitle('Treatment effects within each timepoint')+ 
  ylim(-.1, .3)




# displays estimated means for each treatment group over time.
daily_means <- emmeans(mod_interact, ~ treatment | day) %>% tidy(conf.int=T)


daily_means %>% 
  ggplot(aes(x=day, y=estimate,fill=treatment, ymin=conf.low, ymax=conf.high, group=treatment, color=treatment)) + 
  geom_line(size=1.25, alpha=.9) + 
  geom_errorbar(size=1, width=.2, position=position_dodge(width = .1))+
  geom_point(shape=21, size=4, color='black', position=position_dodge(width = .1)) + 
  theme_cowplot()+ 
  theme(panel.grid.major = element_line(color='grey')) +
  ylab('estimated iRep growth rate')

daily_means %>% 
  ggplot(aes(x=treatment, y=estimate,fill=day, ymin=conf.low, ymax=conf.high, group=day, color=day)) + 
  # geom_line(size=1.25, alpha=.9) + 
  geom_errorbar(size=1, width=.2, position=position_dodge(width = .1))+
  geom_point(shape=21, size=4, color='black', position=position_dodge(width = .1)) + 
  theme_cowplot()+ 
  theme(panel.grid.major = element_line(color='grey')) +
  ylab('estimated iRep growth rate')


### heatmaps?
library(pheatmap)

dat_QC %>% filter(day=='07') %>% 
  transmute(sample, genome,log_rabund=log(`relative abundance`)) %>% 
  pivot_wider(names_from = genome, values_from = log_rabund, values_fill = -1) %>% 
  column_to_rownames(var='sample') %>% 
  pheatmap()
  

color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

dat_QC %>% filter(day=='35') %>% 
  transmute(sample, genome,log_rabund=log(`relative abundance`)) %>% 
  pivot_wider(names_from = genome, values_from = log_rabund, values_fill = -1) %>% 
  column_to_rownames(var='sample') %>% 
  pheatmap()


dat_QC %>% filter(day=='78') %>% 
  transmute(sample, genome,log_rabund=log(`relative abundance`)) %>% 
  pivot_wider(names_from = genome, values_from = log_rabund, values_fill = -1) %>% 
  column_to_rownames(var='sample') %>% 
  pheatmap()


