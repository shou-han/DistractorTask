preloadPackages <- function() {
  if (!require(psych)) {
    install.packages("psych")
  }
  if (!require(nlme)) {
    install.packages("nlme")
  }
  if (!require(car)) {
    install.packages("car")
  }
  if (!require(multcompView)) {
    install.packages("multcompView")
  }
  if (!require(lsmeans)) {
    install.packages("lsmeans")
  }
  if (!require(ggplot2)) {
    install.packages("ggplot2")
  }
  if (!require(rcompanion)) {
    install.packages("rcompanion")
  }
  if (!require(nlme)) {
    install.packages("nlme")
  }
  if (!require(lmerTest)) {
    install.packages("lmerTest")
  }
  if (!require(mediation)) {
    install.packages("mediation")
  }
  if (!require(R.matlab)) {
    install.packages("R.matlab")
  }
  if (!require(TukeyC)) {
    install.packages("TukeyC")
  }
  if (!require(ggpubr)) {
    install.packages("ggpubr")
  }
  ### Install/load required packages
  #List of R packages required for this analysis:
  required_packages <-
    c(
      "psych",
      "ggplot2",
      "tidyr",
      "stringr",
      "lubridate",
      "readxl",
      "knitr",
      "readr",
      "rmarkdown",
      "png",
      "lme4",
      "ez",
      "dplyr",
      "ggpubr"
    )
  #Install required_packages:
  new.packages <-
    required_packages[!(required_packages %in% installed.packages()[, "Package"])]
  if (length(new.packages))
    install.packages(new.packages)
  #Load required_packages:
  lapply(required_packages, require, character.only = TRUE)
  
  
  #Set decimal points and disable scientific notation
  options(digits = 3, scipen = 999)
}
genData <- function(dataRaw) {
  indx = which(!is.na(dataRaw[1, ]))
  a <- dataRaw[1, !is.na(dataRaw[1, ])]
  #omit the first row
  dataRaw = dataRaw[-1, ]
  dataRaw[, ncol(dataRaw)] <- NULL
  dataProc = matrix(nrow = dim(dataRaw)[1], ncol = dim(dataRaw)[2])
  for (l in 1:ncol(a)) {
    dataProc[, l] <- as.numeric(dataRaw[, l])
  }
  dataAll = data.frame((dataProc), stringsAsFactors = FALSE)
  
  #  dataAll=dataRaw
  #   b<-t(rep(a,each=4))
  #  names(dataAll)<-b
  #  names(dataAll)<-make.unique(names(dataAll))
  names(dataAll) <- a
  return(dataAll)
}

nopermANOVA <- function(dataInterest, interestValue) {
  # the rest
  distdata <- dataInterest %>% group_by(ID, distractor) %>%
    select(ID, distractor, interestValue) %>%
    na.omit() %>%
    gather(key_all, dv,-ID,-distractor)
  distMean <- tapply(distdata$dv, distdata$distractor, mean)
  # the binsdata
  bindata <- dataInterest %>% group_by(ID, distractor) %>%
    select(ID, distractor, interestValue) %>%
    na.omit() %>%
    gather(key_all, dv,-ID,-distractor)
  
  #datapleveldist,dataplevelbin,dataplevelrest,
  # the rest
  # resdata<- dataInterest %>% group_by(ID,distractor,bins,hand,Hemi)%>%
  #    select(ID, distractor, interestValue, hand,Hemi, bins) %>%
  #     na.omit() %>%
  #     gather(key_all, dv, -ID, -Hemi, -distractor, -bins,-hand)
  distractData <- ezANOVA(
    data = distdata
    ,
    dv = .(dv)
    ,
    wid = .(ID)
    ,
    within = .(distractor)
    ,
    within_full = .(distractor)
    ,
    within_covariates = NULL
    ,
    between = NULL
    ,
    between_covariates = NULL
    ,
    observed = NULL
    ,
    type = 3
  )
  statData <- ezStats(
    data = distdata
    ,
    dv = .(dv)
    ,
    wid = .(ID)
    ,
    within = .(distractor)
    ,
    within_full = .(distractor)
    ,
    within_covariates = NULL
    ,
    between = NULL
    ,
    between_covariates = NULL
    ,
    type = 3
  )
  # allstatData  <- ezANOVA(data = resdata
  #                         , dv = .(dv)
  #                         , wid = .(ID)
  #                         , within = .(Hemi,distractor,bins,hand)
  #                         , within_covariates = NULL
  #                         , between = NULL
  #                         , between_covariates = NULL
  #                         , observed = NULL
  #                         , type = 3)
  #  binMean <- tapply(distdata$dv, list(bindata$distractor, bindata$bins),mean)
  #  binData<- ezANOVA(data = bindata
  #                          , dv = .(dv)
  #                          , wid = .(ID)
  #                          , within = .(distractor,bins)
  #                          , within_full = .(distractor,bins)
  #                          , within_covariates = NULL
  #                          , between = NULL
  #                          , between_covariates = NULL
  #                          , observed = NULL
  #                          , type = 3)
  return(list(distMean, distractData, statData))#,allstatData))
}

permANOVA <- function(dataInterest, interestValue, permNo) {
  # the rest
  distdata <- dataInterest %>% group_by(ID, distractor) %>%
    select(ID, distractor, interestValue) %>%
    na.omit() %>%
    gather(key_all, dv,-ID,-distractor)
  # the binsdata
  bindata <- dataInterest %>% group_by(ID, distractor) %>%
    select(ID, distractor, interestValue) %>%
    na.omit() %>%
    gather(key_all, dv,-ID,-distractor)
  #datapleveldist,dataplevelbin,dataplevelrest,
  # the rest
  # resdata<- dataInterest %>% group_by(ID,distractor,bins,hand,Hemi)%>%
  #     select(ID, distractor, interestValue, hand,Hemi, bins) %>%
  #     na.omit() %>%
  #     gather(key_all, dv, -ID, -Hemi, -distractor, -bins,-hand)
  #  log <- capture.output({distractData <- ezPerm(data = distdata
  #                         , dv = .(dv)
  #                         , wid = .(ID)
  #                         , within = .(distractor)
  #                         , perms = permNo)});
  #  print('distractDone')
  #  log<- capture.output({  binData <- ezPerm(data = bindata
  #                           , dv = .(dv)
  #                           , wid = .(ID)
  #                           , within = .(distractor,bins)
  #                           , perms = permNo)});
  # print('binDone')
  #  log <- capture.output({  allstatData <- ezPerm(data = resdata
  #                           , dv = .(dv)
  #                           , wid = .(ID)
  #                           , within = .(Hemi,distractor,bins, hand)
  #                           , perms = permNo)});
  # print('alldone')
  return()#list(distractData,binData))#,allstatData))
}

nopermANOVAO <- function(dataInterest, interestValue) {
  # the rest
  distdata <- dataInterest %>% group_by(ID, distractor) %>%
    select(ID, distractor, interestValue) %>%
    na.omit() %>%
    gather(key_all, dv,-ID,-distractor)
  # the binsdata
  bindata <- dataInterest %>% group_by(ID, distractor) %>%
    select(ID, distractor, interestValue) %>%
    na.omit() %>%
    gather(key_all, dv,-ID,-distractor)
  distractData <- ezANOVA(
    data = distdata
    ,
    dv = .(dv)
    ,
    wid = .(ID)
    ,
    within = .(distractor)
    ,
    within_covariates = NULL
    ,
    between = NULL
    ,
    between_covariates = NULL
    ,
    observed = NULL
    ,
    type = 3
  )
  binData <- ezANOVA(
    data = bindata
    ,
    dv = .(dv)
    ,
    wid = .(ID)
    ,
    within = .(distractor)
    ,
    within_covariates = NULL
    ,
    between = NULL
    ,
    between_covariates = NULL
    ,
    observed = NULL
    ,
    type = 3
  )
  return(list(distractData, binData))
}

plotN2figure <-
  function(datapoints,
           interestValue,
           xlabelT,
           ylabelT,
           ylimits) {
    dataPlot <- datapoints %>% group_by(ID, distractor) %>%
      select(ID, distractor, interestValue) %>% na.omit() %>%
      gather(key_all, dv,-ID,-distractor) %>%
      summarize(minInt = mean(-dv, na.rm = TRUE))
    summaryPlot <- dataPlot %>% group_by(distractor) %>%
      summarize(stdInt = se(minInt),
                minInt = mean(minInt, na.rm = TRUE))
    
    plotInt <-
      summaryPlot %>% ggplot() + aes(
        x = distractor,
        y = minInt,
        color = distractor,
        fill = distractor
      ) +
      #geom_jitter(position = position_jitterdodge(jitter.width=0.5,dodge.width=0.9),aes(color=distractor),alpha=0.3)+
      geom_pointrange(
        aes(ymin = minInt - stdInt, ymax = minInt + stdInt),
        position = position_dodge(0.9),
        data = summaryPlot
      ) +
      geom_col(width = 0.5, alpha = 0.3) +
      #  geom_boxplot( width=0.15,fatten=NULL)+
      xlab(xlabelT) + ylab(ylabelT) +
      theme(
        axis.title.x = element_text(face = "bold", size = 20),
        axis.text.x  = element_text(
          face = "bold",
          angle = 0,
          size = 18
        ),
        axis.title.y = element_text(face = "bold", size = 20),
        axis.text.y  = element_text(
          face = "plain",
          angle = 0,
          size = 18
        ),
        legend.text  = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = "none",
        panel.background = element_blank()
      ) +
      scale_x_discrete(labels = c("Absent", "Present")) +
      scale_color_discrete((name = xlabelT), labels = c("Absent", "Present")) +
      coord_cartesian(ylim = ylimits) +
      geom_hline(
        yintercept = 0,
        color = "black",
        size = 0.5,
        alpha = 0.5
      ) +
      scale_fill_manual(values = colorsUsed, name = xlabelT)  +
      scale_color_manual(values = colorsUsed, name = xlabelT)
    # stat_compare_means(data=dataPlot, paired=TRUE, method="t.test",comparisons=list(c("DA","DP")),label="p.signif")
    return(list(plotInt))
  }
# a function for calculating standard error
se=function(x) sqrt(var(x, na.rm=TRUE)/20)


mediationStats <-
  function(dataInterest,
           treatsValue,
           mediateValue,
           outValue) {
    modelMed.O <-
      lm(as.formula(paste(outValue, "~", treatsValue)), dataInterest)
    modelMed.M <-
      lm(as.formula(paste(mediateValue, "~", treatsValue)), dataInterest)
    modelMed.Y <-
      lm(as.formula(paste(
        outValue, "~", mediateValue, "+", treatsValue
      )), dataInterest)
    results <-
      mediate(
        modelMed.M,
        modelMed.Y,
        treat = paste(treatsValue),
        mediator = paste(mediateValue),
        boot = TRUE,
        sims = 500
      )
    return(list(
      modelMed.O,
      modelMed.M,
      modelMed.Y,
      summary(results),
      plot(results)
    ))
    
  }
Normality_tests <- function(data) {
  data %>% ungroup() %>%
    select("ID", "distractor", contains("mean")) %>%
    gather(key, value, -ID, -distractor) %>%
    group_by(key) %>%
    do(
      ShapiroWilk_p_value = shapiro.test(.$value)[2],
      Anderson_Darling_p_value = nortest::ad.test(.$value)[2],
      CramerVonMises_p_value = nortest::cvm.test(.$value)[2],
      Shapiro_Francia_p_value = nortest::sf.test(.$value)[2],
      Kolmogorov_Smirnov_p_value = nortest::lillie.test(.$value)[2]
    ) %>%
    mutate(
      ShapiroWilk_p_value = unlist(ShapiroWilk_p_value),
      Anderson_Darling_p_value = unlist(Anderson_Darling_p_value),
      CramerVonMises_p_value = unlist(CramerVonMises_p_value),
      Shapiro_Francia_p_value = unlist(Shapiro_Francia_p_value),
      Kolmogorov_Smirnov_p_value = unlist(Kolmogorov_Smirnov_p_value)
    ) %>%
    mutate(
      average_p_value = (
        ShapiroWilk_p_value +
          Anderson_Darling_p_value +
          CramerVonMises_p_value +
          Shapiro_Francia_p_value +
          Kolmogorov_Smirnov_p_value
      ) / 5
    ) %>% arrange(average_p_value)
  kable(
    Normality_tests,
    format.args = list(big.mark = ","),
    digits = 4,
    
    caption = "Normality tests"
  )
}
