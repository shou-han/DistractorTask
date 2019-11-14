## this is the stats functions used by distractor data
getSpeedStats<-function(data,test_var,...){
  group_var <- enquos(...)
  summary_var<-enquo(test_var)
  Speed_checker <- data %>% ungroup()%>% group_by(!!! group_var) %>%
    summarise(meanAll = mean(!! summary_var , na.rm=TRUE))
  summary(Speed_checker)
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                  , dv = .(meanAll)
                  , wid = .(ID)
                  , within = .(distractor) 
                  , within_covariates = NULL
                  , between = NULL
                  , between_covariates = NULL
                  , observed = NULL
                  , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);

  model<-lmer(meanAll~distractor+Hemi+distractor*Hemi+(1+distractor+Hemi|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor)
                    , within_covariates = NULL
                    , between = NULL
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}

######################################################################################################################
getAccuracyStats<-function(data,filtervar,test_var,...){
  group_var <- enquos(...)
  summary_var<-enquo(test_var)
  Speed_checker <- data %>% filter(Accuracy==filtervar)%>% ungroup()%>% group_by(!!! group_var) %>%
    summarise(meanAll = mean(!! summary_var , na.rm=TRUE))
  summary(Speed_checker)
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                   , dv = .(meanAll)
                   , wid = .(ID)
                   , within = .(distractor)
                   , within_covariates = NULL
                   , between = NULL
                   , between_covariates = NULL
                   , observed = NULL
                   , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);
  
  model<-lmer(meanAll~distractor+(1|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor)
                    , within_covariates = NULL
                    , between = NULL
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}

######################################################################################################################
getciStats<-function(data,test_var,...){
  group_var <- enquos(...)
  summary_var<-enquo(test_var)
  Speed_checker <- data %>% ungroup()%>% group_by(!!! group_var) %>%
    summarise(meanAll = mean(!! summary_var , na.rm=TRUE))
  summary(Speed_checker)
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                   , dv = .(meanAll)
                   , wid = .(ID)
                   , within = .(distractor,contraipsi) 
                   , within_covariates = NULL
                   , between = NULL
                   , between_covariates = NULL
                   , observed = NULL
                   , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);
  
  model<-lmer(meanAll~distractor+contraipsi+distractor*contraipsi+(1+distractor+contraipsi|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor*contraipsi), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor,contraipsi)
                    , within_covariates = NULL
                    , between = NULL
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}
######################################################################################################################


getExpStats<-function(data,test_var,...){
  group_var <- enquos(...)
  summary_var<-enquo(test_var)
  Speed_checker <- data %>% ungroup()%>% group_by(!!! group_var) %>%
    summarise(meanAll = mean(!! summary_var , na.rm=TRUE))
  summary(Speed_checker)
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                   , dv = .(meanAll)
                   , wid = .(ID)
                   , within = .(distractor) 
                   , within_covariates = NULL
                   , between = .(experiment)
                   , between_covariates = NULL
                   , observed = NULL
                   , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);
  print("########################################################################\n")
  print(paste("Other Statistics Test for", summary_var[2]))
  print("########################################################################\n")  
  model<-lmer(meanAll~experiment+distractor+experiment*distractor+(1|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor)
                    , within_covariates = NULL
                    , between = .(experiment)
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}

######################################################################################################################


getStats<-function(data,test_var){
  summary_var <- enquo(test_var)
  Speed_checker <- data %>% mutate(meanAll = (!! summary_var))
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                   , dv = .(meanAll)
                   , wid = .(ID)
                   , within = .(distractor) 
                   , within_covariates = NULL
                   , between = NULL
                   , between_covariates = NULL
                   , observed = NULL
                   , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);
  
  model<-lmer(meanAll~distractor+(1|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor)
                    , within_covariates = NULL
                    , between = NULL
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}

## this is the stats functions used by distractor data
getSpeedCongStats<-function(data,test_var,...){
  group_var <- enquos(...)
  summary_var<-enquo(test_var)
  Speed_checker <- data %>% ungroup()%>% group_by(!!! group_var) %>%
    summarise(meanAll = mean(!! summary_var , na.rm=TRUE))
  summary(Speed_checker)
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                   , dv = .(meanAll)
                   , wid = .(ID)
                   , within = .(distractor,cong) 
                   , within_covariates = NULL
                   , between = NULL
                   , between_covariates = NULL
                   , observed = NULL
                   , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);
  
  model<-lmer(meanAll~distractor+cong+distractor*cong+(1+distractor+cong|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor)
                    , within_covariates = NULL
                    , between = NULL
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}


######################################################################################################################
getAccuracyCongStats<-function(data,filtervar,test_var,...){
  group_var <- enquos(...)
  summary_var<-enquo(test_var)
  Speed_checker <- data %>% filter(Accuracy==filtervar)%>% ungroup()%>% group_by(!!! group_var) %>%
    summarise(meanAll = mean(!! summary_var , na.rm=TRUE))
  summary(Speed_checker)
  perRT <- ezANOVA(data = data.frame(Speed_checker)
                   , dv = .(meanAll)
                   , wid = .(ID)
                   , within = .(distractor,cong)
                   , within_covariates = NULL
                   , between = NULL
                   , between_covariates = NULL
                   , observed = NULL
                   , type=3);
  print("########################################################################\n")
  print(paste("Statistics Test for", summary_var[2]))
  print("########################################################################\n")
  print(perRT);
  
  model<-lmer(meanAll~distractor+cong+distractor*cong+(1+distractor+cong|ID),data=data.frame(Speed_checker), REML=TRUE)
  print(anova(model))
  print(emmeans(model, list(pairwise ~ distractor), adjust = "bonf"))
  
  
  statRT <- ezStats(data = Speed_checker
                    , dv = .(meanAll)
                    , wid = .(ID)
                    , within = .(distractor,cong)
                    , within_covariates = NULL
                    , between = NULL
                    , between_covariates = NULL
                    , type = 3)
  print("########################################################################\n")
  print(paste("Statistics for", summary_var[2]))
  print("########################################################################\n")
  print(statRT);
}
