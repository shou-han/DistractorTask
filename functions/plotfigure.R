plotfigurebar<-function(data, xvar,yvar,groupvar,xlabelT,ylabelT,legendlabel,colorsUsed,ylims,filtervar,filename){
  yvar <- enquo(yvar)
  xvar<-enquo(xvar)
  groupvar<-enquo(groupvar)
  
  Speed_checker_part <- data %>% group_by(ID,!! xvar, !! groupvar) %>% summarise(meanyvarT=mean(!! yvar, na.rm=TRUE))
 
   Speed_checker_summary<- Speed_checker_part%>% group_by(!! xvar, !! groupvar) %>%
    summarise(meanyvar  = mean(meanyvarT, na.rm=TRUE), stdyvar = se(meanyvarT))
  dataplot<- Speed_checker_summary  %>% filter(Accuracy==filtervar)

  #png(filename,width=400,heigh=350, units="px")
  
  plotSpeed <- dataplot %>% ggplot()+aes(x=!! xvar, y=meanyvar, color=!! xvar, fill=!! xvar) +
    geom_pointrange(aes(ymin=meanyvar-stdyvar,ymax=meanyvar+stdyvar),position = position_dodge(0.9))+
    geom_col(width=0.5,alpha=0.3)+
    xlab(xlabelT) + ylab(ylabelT) +
    theme(axis.title.x = element_text(face="bold", size=25),
          axis.text.x  = element_text(angle=0,  size=25),
          axis.title.y = element_text(face="bold", size=25),
          axis.text.y  = element_text(angle=0, size=25),
          legend.text  =element_text(size=25),
          legend.title = element_text(size=25),
          legend.position="none",
          panel.background = element_blank())+
    scale_x_discrete(labels=legendlabel)+
    scale_color_discrete((name=xlabelT),labels=legendlabel)+
    scale_fill_manual(values=colorsUsed,name=xlabelT, labels=legendlabel)  +
    scale_color_manual(values=colorsUsed,name=xlabelT, labels=c("Absent", "Present"))+
    coord_cartesian(ylim=ylims)#c(800,1000))+  
  #plotSpeed
  #dev.off()
  ggsave(file=filename, dpi=300,width=170,height=150, units="mm")
  print(plotSpeed)
  return(Speed_checker_summary)
}

plotN2bar<-function(data, xvar,yvar,groupvar,xlabelT,ylabelT,legendlabel,colorsUsed,ylims,filename){
  yvar <- enquo(yvar)
  xvar<-enquo(xvar)
  groupvar<-enquo(groupvar)
  
  Speed_checker_part <- data %>% group_by(ID,!! xvar, !! groupvar) %>% summarise(meanyvarT=mean(!! yvar, na.rm=TRUE))
  
  Speed_checker_summary<- Speed_checker_part%>% group_by(!! xvar, !! groupvar) %>%
    summarise(meanyvar  = mean(meanyvarT, na.rm=TRUE), stdyvar = se(meanyvarT))
  dataplot<- Speed_checker_summary
  plotSpeed <- dataplot %>% ggplot()+aes(x=!! xvar, y=meanyvar, color=!! xvar, fill=!! xvar) +
    geom_pointrange(aes(ymin=meanyvar-stdyvar,ymax=meanyvar+stdyvar),position = position_dodge(0.9))+
    geom_col(width=0.5,alpha=0.3)+
    xlab(xlabelT) + ylab(ylabelT) +
    theme(axis.title.x = element_text(face="bold", size=25),
          axis.text.x  = element_text(angle=0,  size=25),
          axis.title.y = element_text(face="bold", size=25),
          axis.text.y  = element_text(angle=0, size=25),
          legend.text  =element_text(size=25),
          legend.title = element_text(size=25),
          legend.position="none",
          panel.background = element_blank())+
    geom_hline(yintercept=0, color = "black", size=0.5,alpha=0.5)+
    scale_x_discrete(labels=legendlabel)+
    scale_color_discrete((name=xlabelT),labels=legendlabel)+
    scale_fill_manual(values=colorsUsed,name=xlabelT, labels=legendlabel)  +
    scale_color_manual(values=colorsUsed,name=xlabelT, labels=c("Absent", "Present"))+
    coord_cartesian(ylim=ylims)#c(800,1000))+  
  ggsave(file=filename, dpi=300, width=170,height=150, units="mm")
  print(plotSpeed)
  
  return(Speed_checker_summary)
}
plotviolin<-function(data,yvar,xvar,groupvar,xlabelT, ylabelT,legendlabel,filename){
  yvar <- enquo(yvar)
  xvar<-enquo(xvar)
  groupvar<-enquo(groupvar)
  ANOVA_dataPlot<-  data%>% 
    group_by(ID, !! xvar, !! groupvar) %>%
    summarise(meanyvar  = mean(!! yvar, na.rm=TRUE), stdyvar = se(!! yvar))
  
  plotSpeed <- ANOVA_dataPlot %>% ggplot()+aes(x=!! xvar, y=meanyvar, color=!! groupvar) +
    geom_violin() +
    geom_boxplot(alpha=0.1) +
    geom_point() +
    # facet_grid(~hand, labeller = "label_both")+  
    geom_point() +
    xlab(xlabelT) + ylab(ylabelT) +
    theme(axis.title.x = element_text(face="bold", size=25),
          axis.text.x  = element_text(face="bold", angle=0,  size=25),
          axis.title.y = element_text(face="bold", size=25),
          axis.text.y  = element_text(face="bold", angle=0, size=25),
          legend.text  =element_text(size=25),
          panel.background = element_blank())+
    labs(color=legendlabel)
  ggsave(file=filename, dpi=300,width=170,height=150, units="mm")
  print(plotSpeed)
}

plotCPPbar<-function(data, xvar,yvar,groupvar,xlabelT,ylabelT,legendlabel,colorsUsed,ylims,filename){
  yvar <- enquo(yvar)
  xvar<-enquo(xvar)
  groupvar<-enquo(groupvar)
  
  Speed_checker_part <- data %>% group_by(ID,!! xvar, !! groupvar) %>% summarise(meanyvarT=mean(!! yvar, na.rm=TRUE))
  
  Speed_checker_summary<- Speed_checker_part%>% group_by(!! xvar, !! groupvar) %>%
    summarise(meanyvar  = mean(meanyvarT, na.rm=TRUE), stdyvar = se(meanyvarT))
  dataplot<- Speed_checker_summary

  plotSpeed <- dataplot %>% ggplot()+aes(x=!! xvar, y=meanyvar, color=!! xvar, fill=!! xvar) +
    geom_pointrange(aes(ymin=meanyvar-stdyvar,ymax=meanyvar+stdyvar),position = position_dodge(0.9))+
    geom_col(width=0.5,alpha=0.3)+
    xlab(xlabelT) + ylab(ylabelT) +
    theme(axis.title.x = element_text(face="bold", size=25),
          axis.text.x  = element_text(angle=0,  size=25),
          axis.title.y = element_text(face="bold", size=25),
          axis.text.y  = element_text(angle=0, size=25),
          legend.text  =element_text(size=25),
          legend.title = element_text(size=25),
          legend.position="none",
          panel.background = element_blank())+
    geom_hline(yintercept=0, color = "black", size=0.5,alpha=0.5)+
    scale_x_discrete(labels=legendlabel)+
    scale_color_discrete((name=xlabelT),labels=legendlabel)+
    scale_fill_manual(values=colorsUsed,name=xlabelT, labels=legendlabel)  +
    scale_color_manual(values=colorsUsed,name=xlabelT, labels=c("Absent", "Present"))+
    scale_y_continuous(expand = c(0, 0),limits=ylims)
  ggsave(file=filename, dpi=300,width=170,height=150, units="mm")
  print(plotSpeed)
  
  return(Speed_checker_summary)
}