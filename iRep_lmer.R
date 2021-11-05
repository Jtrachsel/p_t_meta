
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

LOOK |> ggplot(aes(x=iRep, y=resids)) + geom_point()

treatment_effects <- emmeans(mod, ~ treatment | day) %>%
  contrast(method='pairwise', adjust='fdr')

time_effects <- emmeans(mod, ~ day | treatment) %>%
  contrast(method='pairwise', adjust='fdr') 


all_effects <- rbind(treatment_effects, time_effects) %>%tidy(conf.int=T) |> 
  mutate(contrast=factor(contrast, levels = ''))


library(cowplot)
# treatment effects at different timepoints
all_effects |>
  filter(treatment != '.') |>
  ggplot(aes(x=contrast,color=treatment, y=estimate, ymin=conf.low, ymax=conf.high)) + 
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
  ggplot(aes(x=contrast,color=day, y=estimate, ymin=conf.low, ymax=conf.high)) + 
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
daily_means <- emmeans(mod, ~ treatment | day) %>% tidy(conf.int=T)


daily_means %>% 
  ggplot(aes(x=day, y=estimate,fill=treatment, ymin=conf.low, ymax=conf.high, group=treatment, color=treatment)) + 
  geom_line(size=1.25, alpha=.9) + 
  geom_errorbar(size=1, width=.2, position=position_dodge(width = .1))+
  geom_point(shape=21, size=4, color='black', position=position_dodge(width = .1)) + 
  theme_cowplot()+ 
  theme(panel.grid.major = element_line(color='grey'))

  
  

# emmeans(mod, ~  treatment) %>%
#   contrast(method='revpairwise')
# 
# emmeans(mod, ~  day) %>%
#   contrast(method='revpairwise')
# 
# 
# contrasts_of_interest <- function(mod, adjust='fdr'){
#   treatment_effects <- 
#     emmeans(mod, ~ Treatment | WeeksPostWeaning + Tissue) %>%
#     contrast(method='revpairwise')
#   tissue_effects <- 
#     emmeans(mod, ~ Tissue | WeeksPostWeaning + Treatment) %>%
#     contrast(method='revpairwise') 
#   time_effects <- 
#     emmeans(mod, ~ WeeksPostWeaning | Tissue + Treatment) %>%
#     contrast(method='revpairwise')
#   all_comps <- 
#     rbind(treatment_effects,tissue_effects,time_effects, adjust=adjust) %>% 
#     tidy(conf.int=T)
#   
#   return(all_comps)
# }
# 
# models <- 
#   NG10_flowdat %>% 
#   select(AnimalID, WeeksPostWeaning, Tissue,
#          Treatment, contains('pct_totalCD3E')) %>% 
#   pivot_longer(cols =contains('pct_totalCD3E'),
#                names_to = 'cell_type', 
#                values_to = 'pct_tot_CD3E') %>% 
#   group_by(cell_type) %>% 
#   nest() %>% 
#   mutate(
#     simple_mods=
#       map(.x=data, 
#           .f=~lmer(data=.x, 
#                    formula= pct_tot_CD3E ~ Tissue + Treatment + factor(WeeksPostWeaning) + (1|AnimalID))),
#     timeXtissue_interact_mods=
#       map(.x=data, 
#           .f=~lmer(data=.x, 
#                    formula= pct_tot_CD3E ~ Tissue * factor(WeeksPostWeaning) + Treatment + (1|AnimalID))),
#     tissueXtreatment_interact_mods=
#       map(.x=data, 
#           .f=~lmer(data=.x, 
#                    formula= pct_tot_CD3E ~ Tissue *Treatment + factor(WeeksPostWeaning)  + (1|AnimalID))),
#     timeXtreatment_interact_mods=
#       map(.x=data, 
#           .f=~lmer(data=.x, 
#                    formula= pct_tot_CD3E ~ Treatment * factor(WeeksPostWeaning) + Tissue + (1|AnimalID))),
#     mods=
#       map(.x=data, 
#           .f=~lmer(data=.x,
#                    formula = pct_tot_CD3E ~ Tissue * Treatment * factor(WeeksPostWeaning) + (1|AnimalID))))
# 
# anova(models[1,3][[1]][[1]],models[1,4][[1]][[1]], models[1,5][[1]][[1]], models[1,6][[1]][[1]], models[1,7][[1]][[1]])
# 
# 
# models$simple_mods %>% map(summary)
# models$timeXtissue_interact_mods %>% map(summary)
# models$tissueXtreatment_interact_mods%>% map(summary)
# models$timeXtreatment_interact_mods%>% map(summary)
# models$mods%>% map(summary)
# 
# models_and_contrasts$simple_mods[[1]] %>% summary()
# models_and_contrasts$simple_mods[[2]] %>% summary()
# models_and_contrasts$simple_mods[[3]] %>% summary()
# models_and_contrasts$simple_mods[[4]] %>% summary()
# 
# 
# models_and_contrasts$time_interact_mods[[1]] %>% summary()
# models_and_contrasts$time_interact_mods[[2]] %>% summary()
# models_and_contrasts$time_interact_mods[[3]] %>% summary()
# models_and_contrasts$time_interact_mods[[4]] %>% summary()
# 
# 
# models_and_contrasts$mods[[1]] %>% summary()
# models_and_contrasts$mods[[2]] %>% summary()
# models_and_contrasts$mods[[3]] %>% summary()
# models_and_contrasts$mods[[4]] %>% summary()
# 
# 
# models_and_contrasts <- 
#   models %>% 
#   mutate(
#     comps=map(.x=mods, .f=~contrasts_of_interest(.x)), 
#     time_interact_comps=map(.x=time_interact_mods, .f=~contrasts_of_interest(.x)),
#     simple_model_comps=map(.x=simple_mods, .f=~contrasts_of_interest(.x))) 
# 
# 
# models_and_contrasts$simple_mods[[1]] %>% summary()
# models_and_contrasts$simple_mods[[2]] %>% summary()
# models_and_contrasts$simple_mods[[3]] %>% summary()
# models_and_contrasts$simple_mods[[4]] %>% summary()
# 
# 
# models_and_contrasts$time_interact_mods[[1]] %>% summary()
# models_and_contrasts$time_interact_mods[[2]] %>% summary()
# models_and_contrasts$time_interact_mods[[3]] %>% summary()
# models_and_contrasts$time_interact_mods[[4]] %>% summary()
# 
# 
# models_and_contrasts$mods[[1]] %>% summary()
# models_and_contrasts$mods[[2]] %>% summary()
# models_and_contrasts$mods[[3]] %>% summary()
# models_and_contrasts$mods[[4]] %>% summary()
# 
# 
# models_and_contrasts$simple_model_comps
# 
# # Model with all possible interactions
# 
# 
# all_contrasts <- 
#   models_and_contrasts %>% 
#   select(cell_type, comps) %>% 
#   unnest(cols = comps) %>% 
#   mutate(effect_type=
#            case_when(
#              grepl('um', contrast) ~ 'Tissue',
#              grepl('Carb', contrast) ~ 'Abx', 
#              grepl('4', contrast)    ~ 'Time'
#            ), 
#          `Cell Type` = sub('_pct_totalCD3E','',cell_type))
# 
# ### simple model ###
# 
# simple_contrasts <- 
#   models_and_contrasts %>% 
#   select(cell_type, simple_model_comps) %>% 
#   unnest(cols = simple_model_comps) %>% 
#   mutate(effect_type=
#            case_when(
#              grepl('um', contrast) ~ 'Tissue',
#              grepl('Carb', contrast) ~ 'Abx', 
#              grepl('4', contrast)    ~ 'Time'
#            ), 
#          `Cell Type` = sub('_pct_totalCD3E','',cell_type))
# 
# # Time interaction
# 
# time_inter_contrasts <- 
#   models_and_contrasts %>% 
#   select(cell_type, time_interact_comps) %>% 
#   unnest(cols = time_interact_comps) %>% 
#   mutate(effect_type=
#            case_when(
#              grepl('um', contrast) ~ 'Tissue',
#              grepl('Carb', contrast) ~ 'Abx', 
#              grepl('4', contrast)    ~ 'Time'
#            ), 
#          `Cell Type` = sub('_pct_totalCD3E','',cell_type))
# 
# 
# ### By effect types
# 
# # All interactions
# 
# all_contrasts %>%
#   filter(effect_type == 'Abx') %>% 
#   ggplot(aes(x=Tissue, y=estimate, ymin=conf.low, ymax=conf.high, color=WeeksPostWeaning)) + 
#   geom_pointrange(position = position_dodge(width = .25)) + geom_hline(yintercept = 0)+
#   facet_wrap(~`Cell Type`, nrow = 1, scales = 'free_x') + 
#   coord_flip() + 
#   ylab('Difference between Carbadox and NoAbx (Carb-NoAbx, % tot CD3E)') + 
#   ggtitle('Abx Effects') +
#   theme_bw()
# 
# 
# all_contrasts %>%
#   filter(effect_type == 'Tissue') %>% 
#   ggplot(aes(x=contrast, y=estimate, ymin=conf.low, ymax=conf.high, color=WeeksPostWeaning, shape=Treatment)) + 
#   geom_pointrange(position = position_dodge(width = .35)) + geom_hline(yintercept = 0)+
#   facet_wrap(~`Cell Type`, nrow = 1, scales = 'free_x') + 
#   coord_flip() + 
#   ylab('Difference between Tissues (% tot CD3E)') + 
#   ggtitle('Tissue Effects')+
#   theme_bw()
