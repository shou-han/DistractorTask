## this has preloading of Nicole's data
preloadData <- function(dataStats) {
  dataF <- dataStats
  dataF$log_RT <- log(dataF$RT) #log
  #####Z-score each participant's log_RT data ####
  dataF$IDbyITIbyHemifieldbygroup <-
    interaction(dataF$Subject, dataF$hemi, dataF$group)
  m <-
    tapply(dataF$log_RT, dataF$IDbyITIbyHemifieldbygroup, mean, na.rm = T)
  s <-
    tapply(dataF$log_RT, dataF$IDbyITIbyHemifieldbygroup, sd, na.rm = T)
  #calculate log_RT.Z and save it inside data.frame
  dataF$logz1 <- (dataF$log_RT - m[dataF$IDbyITIbyHemifieldbygroup])
  dataF$logz1 <- s[dataF$IDbyITIbyHemifieldbygroup]
  dataF$log_RT.Z <-
    (dataF$log_RT - m[dataF$IDbyITIbyHemifieldbygroup]) / s[dataF$IDbyITIbyHemifieldbygroup]
  dataF <-
    dataF %>% dplyr::rename(
      ID = Subject,
      Hemi = hemi,
      CPPlevel = CPPampatresponse,
      N2cLatency = N2c_T,
      N2iLatency = N2i_T,
      N2cpeak = N2cP,
      N2ipeak = N2iP,
      CPPpeakTime = CPP_T
    )
  dataF <- dataF %>% #make factors
    mutate(
      Hemi = ifelse(Hemi == 1, "Left", "Right"),
      hand = ifelse(hand == 1, "Left", "Right")
    ) %>%
    mutate_each_(funs(factor), c("group", "Hemi", "ID")) #make factor
  levels(dataF$group) <- list("AgeM" = "1",
                              "ReadM" = "2",
                              "Dyslexic" = "3")
  #Remove trials where absolute log_RT.Z>3 (i.e. remove outlier RTs)
  #dataF<-dataF[!abs(dataF$log_RT.Z)>5,]
  return(dataF)
}

preloadSpeedData <- function(dataA) {
  data <- dataA %>% #make factors
    mutate(
      side = ifelse(side == 1, "Left", "Right"),
      hand = ifelse(hand == 1, "left", "right"),
      hit = ifelse(is.na(hit), "miss", ifelse(hit == 1, "hit", "wrong"))
    ) %>% #rename
    rename(
      .,
      ID = subjectAll,
      Accuracy = hit,
      Hemi = side,
      distractor = "c"
    ) %>%
    mutate_each_(funs(factor), c("distractor", "Hemi", "iti", "ID", "hand")) #make factor
  levels(data$distractor) <- list("DA" = "1", "DP" = "2")
  data$RT = data$RT * 2
  data <- filter(data, RT < 1500, RT > 150)
  
  data$log_RT <- (data$RT) #log
  #####Z-score each participant's log_RT data ####
  data$IDbyITIbyHemifield <- interaction(data$ID, data$iti, data$hand)
  #calculate mean and sd
  m <- tapply(data$log_RT, data$IDbyITIbyHemifield, mean, na.rm = T)
  s <- tapply(data$log_RT, data$IDbyITIbyHemifield, sd, na.rm = T)
  #calculate log_RT.Z and save it inside data.frame
  data$log_RT.Z <-
    as.numeric((data$log_RT - m[data$IDbyITIbyHemifield]) / s[data$IDbyITIbyHemifield])
  #Remove trials where absolute log_RT.Z>3 (i.e. remove outlier RTs)
  data <- data[!abs(data$log_RT.Z) > 3, ]
  return(data)
}

preloadAccuracydata <- function(dataA) {
  data <- dataA %>% #make factors
    mutate(
      side = ifelse(side == 1, "Left", "Right"),
      c = ifelse(c == 2, "DP", "DA"),
      hand = ifelse(hand == 1, "left", "right"),
      hit = ifelse(is.na(hit), "miss", ifelse(hit == 1, "hit", "wrong"))
    ) %>% #rename
    rename(
      .,
      ID = subjectAll,
      Accuracy = hit,
      Hemi = side,
      distractor = c
    ) %>%
    mutate_each_(funs(factor), c("distractor", "Hemi", "iti", "ID", "hand")) #make factor
  data <- data %>% group_by(ID, distractor) %>%
    mutate(
      hitwronmiss = ifelse(Accuracy == "hit", 1, (ifelse(
        Accuracy == "miss", 0, 1
      ))),
      totalTrials = sum(hitwronmiss)
    ) %>% filter(Accuracy != "miss")
  AccID <- data %>% group_by(ID, distractor, Accuracy, totalTrials) %>%
    summarise(totalNo  = sum(hitwronmiss, na.rm = TRUE)) %>%
    mutate(pect = as.numeric(totalNo / totalTrials * 100))
  return(AccID)
}


preloadAccuracyCongdata <- function(dataA) {
  data <- dataA %>% #make factors
    mutate(
      side = ifelse(side == 1, "Left", "Right"),
      c = ifelse(c == 2, "DP", "DA"),
      hand = ifelse(hand == 1, "left", "right"),
      hit = ifelse(is.na(hit), "miss", ifelse(hit == 1, "hit", "wrong"))
    ) %>% #rename
    rename(
      .,
      ID = subjectAll,
      Accuracy = hit,
      Hemi = side,
      distractor = c
    ) %>%
    mutate_each_(funs(factor), c("distractor", "cong", "iti", "ID", "hand")) #make factor
  data <- data %>% group_by(ID, distractor) %>%
    mutate(
      hitwronmiss = ifelse(Accuracy == "hit", 1, (ifelse(
        Accuracy == "miss", 0, 1
      ))),
      totalTrials = sum(hitwronmiss)
    ) %>% filter(Accuracy != "miss")
  AccID <- data %>% group_by(ID, distractor, Accuracy, totalTrials,cong) %>%
    summarise(totalNo  = sum(hitwronmiss, na.rm = TRUE)) %>%
    mutate(pect = as.numeric(totalNo / totalTrials * 100))
  return(AccID)
}

preloadmissData <- function(dataA) {
  data <- dataA %>% #make factors
    mutate(
      side = ifelse(side == 1, "Left", "Right"),
      c = ifelse(c == 2, "DP", "DA"),
      hand = ifelse(hand == 1, "left", "right"),
      hit = ifelse(is.na(hit), "miss", "hit")
    ) %>% #rename
    rename(
      .,
      ID = subjectAll,
      Accuracy = hit,
      Hemi = side,
      distractor = c
    ) %>%
    mutate_each_(funs(factor), c("distractor", "cong", "iti", "ID", "hand")) #make factor
  testdata <- data %>% group_by(ID, distractor) %>%
    mutate(
      hitwronmiss = ifelse(Accuracy == "hit", 1, 1),
      totalTrials = sum(hitwronmiss)
    )
  #%>%
  #  filter(Accuracy=="hit")
  #perfectSubj = c(1,6 ,7 ,15,13, 17, 18, 20)
  #testdata<- testdata %>% filter(!(ID %in% perfectSubj))# %>% filter(!(Accuracy=="miss"))
  AccID <-
    testdata %>% group_by(ID, distractor, Accuracy, totalTrials) %>%
    summarise(totalNo  = sum(hitwronmiss, na.rm = TRUE)) %>%
    mutate(pect = as.numeric(totalNo / totalTrials * 100)) %>% filter(Accuracy ==
                                                                        "miss")
  return(AccID)
}

preloadN2peakData<- function(plevelAll){
  # Peaks
  ANOVA_data <- plevelAll%>% group_by(ID,distractor)%>%
    select(ID,distractor,meanN2ipeak, meanN2cpeak) %>%
    na.omit() %>%
    gather(key, dv, -ID,-distractor) %>%
    mutate(contraipsi = ifelse(key=="meanN2cpeak","N2c","N2i"))%>%
    mutate_each_(funs(factor),c("contraipsi")) #make factor
  return(ANOVA_data)
}
preloadN2latencyData<- function(plevelAll){
  # Peaks
  ANOVA_data <- plevelAll%>% group_by(ID,distractor)%>%
    select(ID,distractor,meanN2iLatency, meanN2cLatency) %>%
    na.omit() %>%
    gather(key, dv, -ID,-distractor) %>%
    mutate(contraipsi = ifelse(key=="meanN2cLatency","N2c","N2i"))%>%
    mutate_each_(funs(factor),c("contraipsi")) #make factor
  return(ANOVA_data)
}
preloadmissCongData <- function(dataA) {
  data <- dataA %>% #make factors
    mutate(
      side = ifelse(side == 1, "Left", "Right"),
      c = ifelse(c == 2, "DP", "DA"),
      hand = ifelse(hand == 1, "left", "right"),
      hit = ifelse(is.na(hit), "miss", "hit")
    ) %>% #rename
    rename(
      .,
      ID = subjectAll,
      Accuracy = hit,
      Hemi = side,
      distractor = c
    ) %>%
    mutate_each_(funs(factor), c("distractor", "cong", "iti", "ID", "hand")) #make factor
  testdata <- data %>% group_by(ID, distractor) %>%
    mutate(
      hitwronmiss = ifelse(Accuracy == "hit", 1, 1),
      totalTrials = sum(hitwronmiss)
    )
  #%>%
  #  filter(Accuracy=="hit")
  #perfectSubj = c(1,6 ,7 ,15,13, 17, 18, 20)
  #testdata<- testdata %>% filter(!(ID %in% perfectSubj))# %>% filter(!(Accuracy=="miss"))
  AccID <-
    testdata %>% group_by(ID, distractor, Accuracy, cong,totalTrials) %>%
    summarise(totalNo  = sum(hitwronmiss, na.rm = TRUE)) %>%
    mutate(pect = as.numeric(totalNo / totalTrials * 100)) %>% filter(Accuracy ==
                                                                        "miss")
  return(AccID)
}