clear all
close all
set(0,'DefaultFigureVisible','on')
t=1;
addpath(genpath('/scratch/yn70/ShouHan/Distractors/CSDtoolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/eeglab13_6_5b/'));
addpath('function_programs/');
eeglab
namestring ={'150dots','10dots'};
dLabels = {'Absent','Present'};
addpath('function_programs/');
onew=2;
oldnew={'', 'old'};
if onew==2
    ssj =[1:18 20:21];
else
    ssj=[1:20];
end
chn = 13;
CSD = 1;
if CSD
    yCPP = [-10,35]; ytickCPP =[-10:5:35];
    yN2pc = [-30,20]; ytickN2pc = [-30:5:20];
    yN2c= [-30,20]; ytickN2c=[-30:5:20];
    yN2i=[-30,20];ytickN2i=[-30:5:20];
    ylim = [-30 50]; ylimBeta=[-2 0.5];ytickBeta=[-2:0.3:0.5];
else
    yCPP = [-4,8]; ytickCPP =[-4:2:8];
    yN2pc = [-4,4]; ytickN2pc = [-4:2:4];
    yN2c= [-4,4]; ytickN2c=[-4:2:4];
    yN2i=[-4,4];ytickN2i=[-4:2:4];
    ylim = [-4 8];ylimBeta=[-0.1 0.1];ytickBeta=[-0.9:1:3];
end
cs = [1 ;2];
csdiff= [1; 2];
conds=[];
ch_beta = [8 25 8 25];
xlabels_plot = {'\fontsize{12}-100   ', 0', 200:200:1200};
%combine them
for s = 1:size(ssj,2)
    sT = ssj(s);
    if CSD
        load(['Data/ERPs/group_plots_erp_diff_CSD_1binsStats_uB_35Hz_allTrials_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','t','tr','allBins','no_of_bins',...
            'beta_side','beta_r_side','beta_contra','betar_contra','STFT_time','STFT_timer','RT_side',...
            'betaslope','alpha_side','alpha_asym_side','alpha_t','alpha_tr')
    else
        load(['Data/ERPs/group_plots_erp_diff__1binsStats_cong_35Hz_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','t','tr','allBins','no_of_bins',...
            'beta_side','beta_r_side','beta_contra','betar_contra','STFT_time','STFT_timer','RT_side',...
            'betaslope','alpha_side','alpha_asym_side','alpha_t','alpha_tr')
    end
    
    allBins = allBins;
    for bin=1:no_of_bins
        for d=1:length(cs)
            for cong=1
                c = cs(d,:);
                allstuff_a.RTs{s,d,bin} = [allBins.Rts{1,c,:,bin}];
                sizeRTS(s,d) = length([allstuff_a.RTs{s,d,bin}]);
                CPP_sides(s,d,:,:)= squeeze(mean(allBins.CPP_side(:,c,:,:),2));
                CPPr_sides(s,d,:,:) = squeeze(mean(allBins.CPPr_side(:,c,:,:),2));
                N2c_sides(s,d,:,:)=squeeze(mean(allBins.N2c_side(:,c,:,:),2));
                N2i_sides(s,d,:,:) =squeeze(mean(allBins.N2i_side(:,c,:,:),2));
                N2pc_sides(s,d,:,:)=squeeze(mean(allBins.N2pc_side(:,c,:,:),2));
                
                CPPslope_sides(s,d,:,:) =squeeze(mean(allBins.CPPslopes_side(:,c,:,:),2));
                CPPrslope_sides(s,d,:,:) =squeeze(mean(allBins.CPPrslopes_side(:,c,:,:),2));
                
                ERP_sides(s,d,:,:,:) = squeeze(mean(allBins.ERP_side(:,c,:,:,:),2));
                ERPr_sides(s,d,:,:,:) = squeeze(mean(allBins.ERPr_side(:,c,:,:,:),2));
                beta_contra_all(s,d,:,:,:,:) = squeeze(mean(beta_contra(:,c,:,:,:,:),2));
                betar_contra_all(s,d,:,:,:) = squeeze(mean(betar_contra(:,c,:,:,:),2));
                
                beta_side_all(s,d,:,:,:,:) = squeeze(mean(beta_side(:,c,:,:,:,:),2));
                betar_side_all(s,d,:,:,:) = squeeze(mean(beta_r_side(:,c,:,:,:),2));
                
                beta_c_slope_all(s,d,:,:,:,:) = squeeze(mean(betaslope.contra(:,c,:,:,:,:),2));
                beta_i_slope_all(s,d,:,:,:,:) = squeeze(mean(betaslope.ipsi(:,c,:,:,:,:),2));
                betar_c_slope_all(s,d,:,:,:,:) = squeeze(mean(betaslope.contraR(:,c,:,:,:,:),2));
                betar_i_slope_all(s,d,:,:,:,:) = squeeze(mean(betaslope.ipsiR(:,c,:,:,:,:),2));
                
                beta_click_contra(s,d) = mean(mean(beta_contra(:,c,find(STFT_timer==0),:),4),2);
                RT_sides{s,d}= [RT_side{1,c,:}];
                RT_mean(s,d) = mean([RT_sides{s,d}]);
                
                alpha_sides(s,d,:,:) = squeeze(mean(mean(alpha_side(:,c,:,:,:),2),5));
                alpha_asym_sides(s,d,:,:) = squeeze(mean(mean(alpha_asym_side(:,c,:,:,:),2),5));
                
                % ERP_sidesO(s,:,:,:,:,:) = squeeze(allBins.ERP_sideO);
                % ERPr_sidesO(s,:,:,:,:,:) = squeeze(allBins.ERPr_sideO);
            end
        end
    end
end



%% plot reaction time lines
for bin=1:1
    for sj = 1:size(allstuff_a.RTs,1)
        
        for c = 1:2
            Rt_sj_times(sj,c,:) = mean([allstuff_a.RTs{sj,c,:}]);
        end
    end
    RT_group_times(:,:,:) = Rt_sj_times(:,:,:)';
end

%% Colors

cyan        = [0.2 0.8 0.8];

mangenta    = [1 0 1];
lightgreen  = [0 1 0];
red         = rgb('Red');
blue        = [0    0    1.0000];
green       = [0    1.0000    0.4961];
black        = [0 0 0 ];
brown       = rgb('Chocolate');
orange      = rgb('Gold');
purple      = rgb('Magenta');
dred         = [ 0.6953    0.1328    0.1328];
dblue        = [0         0    0.5430];
dgreen       = [0    0.3906         0];
dblack       = [0 0 0];
dbrown       = rgb('Maroon');
dorange      = rgb('OrangeRed');
dpurple      =[0.5 0 1];
Cols = [dgreen;dred;dblue;green;red;blue;dpurple;purple];Cols1=Cols; Cols4=Cols;
ColsBeta=[dbrown;dorange;brown;orange];
ColsRT=[red;blue];1
RTCols=[dblack];
alphas =[0.4 1 0.2 0.6];
alphaRT=[0.4 1];
return
%% Topology Plots
%% ERP scalp topo collapsed across all trials
close all
if CSD
    maplim = [-15 12];
    mapRlim = [-20 20];
    mapDlim = [-10 10];  
else
    maplim = [-1.5 1.5];
    mapRlim = [-2 2];
    mapDlim = [-1 1];
end
%%%%%%%%%%%%%%%%%DA

left_hemi = [1 33 34 4 3 37 36 5 38 6 39 7 9 41 8 40 10 42 11 43 12 15 45 14 44 46 47 16];
right_hemi = [32 63 62 31 30 60 61 27 59 28 58 29 26 56 25 57 21 55 22 54 23 20 51 19 52 50 49 18];
centre_chans = [35 48 2 13 17 64 24 53];

plot_mean = zeros([64,1]);

for c=1
    %N2pc
    hf1 = figure(1)
    subplot(1,1,1)
    plot_mean(left_hemi) = squeeze(mean(mean(mean(mean(ERP_sides(:,c,left_hemi,find(t>=200 & t<400),:)-...
        ERP_sides(:,1,right_hemi,find(t>=200 & t<400),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',left_hemi);
    %N2c/N2i
    hf2=figure(2)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,c,:,find(t>=200 & t<400),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',1:64);
    %N2c/N2i
    hf3=figure(3)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,c,:,find(t>=200 & t<400),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [-15 12], ...
        'electrodes','off','plotchans',1:64);
    %CPP
    hf4=figure(4)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERPr_sides(:,c,:,find(tr>=-100 & tr<100),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
     [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
    %BetaWeird
    hf5=figure (5)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(beta_side_all(:,c,:,find(STFT_time>=600 & STFT_time<800),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [0 6],...%[min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',1:64);
    %Beta
    hf6=figure (6)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(betar_side_all(:,c,:,find(STFT_timer>=-100 & STFT_timer<100),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [0 6], ...
        'electrodes','off','plotchans',1:64);
    %Alpha
    hf7=figure (7)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(alpha_sides(:,c,:,find(alpha_t>=-700 & alpha_t<-50),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',1:64);    
    %Alpha asym
    hf8=figure (8)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(alpha_asym_side(:,c,:,find(alpha_t>=-700 & alpha_t<-50),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',1:64);    
end
saveas(hf1,['Figures/N2pc' namestring{onew} '.png'])
saveas(hf2,['Figures/N2c_i' namestring{onew} '.png'])
saveas(hf3,['Figures/N2i_c' namestring{onew} '.png'])
saveas(hf4,['Figures/CPP' namestring{onew} '.png'])
saveas(hf5,['Figures/Beta' namestring{onew} '.png'])
saveas(hf6,['Figures/betaR' namestring{onew} '.png'])
saveas(hf7,['Figures/alphapre' namestring{onew} '.png'])

%% RT Behaviour
close all
clear hs hb
    for c = 1:2
        temp  =[];
        for s = 1:length(ssj)
            temp = [temp [allstuff_a.RTs{s,c,:}]];
        end
        meanRTs{c} = temp;
        clear temp 
    end
hf = figure(9998)
hold on
cs=0;
for c = [1 2]
    cs=cs+1;
    clear de_temp
    de_temp = [allstuff_a.RTs{:,c,:}];
      set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[400,1600],'xtick',[400:200:1600],'ylim',[0,500],'ytick',[0:150:450],'TickDir','out');%,'ylim',[-1.5,0.5]);
    h = histc(de_temp, [min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)]);
    hs(cs) =plot([min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)],h,'Color',ColsRT(c,:),'Linewidth',1.5);
    hb(cs) = line([mean(meanRTs{c}) mean(meanRTs{c})],[0 500],'Color',ColsRT(c,:),'Linewidth',1.5,'LineStyle','--');
    
end
ylabel({'Frequency'},'FontSize',16,'fontweight','bold')
xlabel('Reaction Time (ms)','FontSize',16,'fontweight','bold')
legend(hs,dLabels,'FontSize',12,'Location','NorthEast','fontweight','bold')
legend boxoff
saveas(hf,['Figures\RT_all_bins' namestring{onew} '.png']);
% make bar graphs
for c=[1 2]
    mean_RTs(c) = mean([allstuff_a.RTs{:,c,:}]);
    std_RTs(c) = std([allstuff_a.RTs{:,c,:}])/sqrt(20);
end
hf = figure(9999)
hold on
cs=0;
for c=1:2
    hs1(c) = bar(c,mean_RTs(c),'Linewidth',1.5);
    set(hs1(c),'FaceColor',ColsRT(c,:));
    hs2(c) = errorbar(c,mean_RTs(c),std_RTs(c)'.');
    hs2(c).Color = dblack;
    hs2(c).LineWidth = 2;
end
ylabel({'Reaction Time (ms)'},'FontSize',16,'fontweight','bold')
axis([0 3 300 1200]);
xticks([1 2]);xticklabels({'Absent','Present'});
set(gca,'Fontsize',20)
yticks(300:300:1200);
legend(hs,dLabels,'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\RT_bar_graphs' namestring{onew} '.png']);

%% N2 graphs
%shadedErrorBar
%N2pc
close all
clear hf hss1
hf = figure(1)
for bin=1:1
    for c=1:2
        N2pc_c = squeeze(mean(mean(N2pc_sides(1:20,c,:,:),4),2));%[2 3 5 6 8 9 11:13 16 17 18 20 1 2 4 7 10 14:15 19 ]
        meanN2pc = mean(N2pc_c,1);
        stdN2pc = std(N2pc_c,1)/sqrt(size(N2pc_c,1));
        hold on
        (c-1)*2+bin
        hss1((bin-1)*2+c) = shadedErrorBar(t,meanN2pc,stdN2pc,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'linewidth',1.5,'FontSize',16,'xlim',[-100,1400],'xTick',[-100,0:200:1200],'XTickLabel',xlabels_plot,'ylim',yN2pc,'ytick',ytickN2pc,'TickDir','out');%,'ylim',[-1.5,0.5]);
        %text(-170,-33.5,'-100','FontSize',12)
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            ylim,'Color',ColsRT(c,:),'LineWidth',1.5,'LineStyle','--');
        ylabel('N2pc (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16','fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine, hss1(2).mainLine],dLabels,'FontSize',12,'Location','SouthEast')
legend boxoff
saveas(hf,['Figures\N2pcxdist' namestring{onew} '.png']);
% N2c 
clear hss3
hf = figure(3)
for bin = 1:1
    hold on
    for c=1:2
        N2c_c_grp = squeeze(mean(mean(ERP_sides(1:20,c,45,:,:,:),2),5));
        meanN2c = mean(N2c_c_grp,1);
        stdN2c = std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[-100,1400],'xtick',[-100, 0:200:1200],'XTickLabel',xlabels_plot,...
            'ylim',yN2c,'ytick',ytickN2c,'TickDir','out');%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2c,'Color',ColsRT(c,:),'LineWidth',1.5,'LineStyle','--');
        ylabel('N2c (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],dLabels,'FontSize',12,'Location','SouthEast')
legend boxoff
saveas(hf,['Figures\N2cxdist' namestring{onew} '.png']);
% N2i 
clear hss3
hf = figure(2)
for bin = 1:1
    hold on
    for c=1:2
        N2i_c_grp = squeeze(mean(mean(ERP_sides(1:20,c,51,:,:,:),2),5));
        meanN2i = mean(N2i_c_grp,1);
        stdN2i = std(N2i_c_grp,1)/sqrt(size(N2i_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],'XTickLabel',xlabels_plot,...
            'ylim',yN2i,'ytick',ytickN2i,'TickDir','out');%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2i,'Color',ColsRT(c,:),'LineWidth',1.5,'LineStyle','--');
        ylabel('N2i (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],dLabels,'FontSize',12,'Location','SouthEast')
legend boxoff
saveas(hf,['Figures\N2ixdist' namestring{onew} '.png']);
%% CPP_CPPr 
close all
clear hss3 hss4
hf = figure(4)
for bin = 1:1
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPP_sides(:,c,:,:),4),2));
        meanCPP = mean(CPP_c,1);
        stdCPP = std(CPP_c,1)/sqrt(size(CPP_c,1));
        clear CPP_c
       % beta_sides = squeeze(mean(mean(beta_side_all(:,c,:),1),5));
        hold on
        (c-1)*2+bin
        hss1((bin-1)*2+c) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1200],'XTickLabel',xlabels_plot,'ylim',yCPP,'ytick',ytickCPP,'TickDir','out');%,'ylim',[-1.5,0.5]);
        ylabel('CPP (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],yCPP,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss3((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylim,...
           'LineWidth',1.5,'LineStyle','--','Color',ColsRT(c,:));
        hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylim,...
           'LineWidth',1.5,'LineStyle','--','Color',ColsRT(c,:));        
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine hss1(2).mainLine],dLabels,'FontSize',12,'Location','SouthEast')
legend boxoff
saveas(hf,['Figures\CPPxdist' namestring{onew} '.png']);
clear hss11
hf = figure(5)
for bin = 1:1
    hold on
    for c=1:2
        CPPr_c = squeeze(mean(mean(CPPr_sides(:,c,:,:),4),2));
        meanCPPr = mean(CPPr_c,1);
        stdCPPr = std(CPPr_c,1)/sqrt(size(CPPr_c,1));
        clear CPPr_c        
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',yCPP,'ytick',ytickCPP,'TickDir','out');%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('CPP (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from Response (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],dLabels,'FontSize',12,'Location','NorthWest')
legend boxoff
saveas(hf,['Figures\CPPrxdist' namestring{onew} '.png']);
%% Pre target Alpha power
alpha_ch = [46 47 48 49 50];
clear hss11
hf = figure(10)
for bin = 1:1
    hold on
    for c=1:2
        alpha_c = squeeze(mean(mean(alpha_sides(:,c,alpha_ch,:),3),2));
        meanalpha = mean(alpha_c,1);
        stdalpha = std(alpha_c,1)/sqrt(size(alpha_c,1));
        clear alpha_c        
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(alpha_t,meanalpha,stdalpha,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0 16],'ytick',[0:2:16],'TickDir','out');%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',ColsRT(c,:),'LineWidth',1.5,'LineStyle','-');
        ylabel('Alpha (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],dLabels,'FontSize',12,'Location','SouthWest')
legend boxoff
saveas(hf,['Figures\alpha_sides' namestring{onew} '.png']);
clear hss11
hf = figure(12)
for bin = 1:1
    hold on
    for c=1:2
        alpha_c = squeeze(mean(mean(alpha_asym_sides(:,c,[19],:),3),2));
        meanalpha = mean(alpha_c,1);
        stdalpha = std(alpha_c,1)/sqrt(size(alpha_c,1));
        clear alpha_c        
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(alpha_t,meanalpha,stdalpha,'lineprops',{'Color',ColsRT(c,:),'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0.1 0.1],'ytick',ytickCPP,'TickDir','out');%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-0.5 0.5],'Color',ColsRT(c,:),'LineWidth',1.5,'LineStyle','-');
        ylabel('Alpha (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],dLabels,'FontSize',12,'Location','NorthWest')
legend boxoff
saveas(hf,['Figures\alpha_sides_asym' namestring{onew} '.png']);
%% beta
clear h1 h2 hss4 hf hss3
% contra beta and ipsi beta x distractors resp
%150dots
%ylimBeta= [2 3];ytickBeta = [2:0.2:3]; %comment out for baselined
%ylimBeta= [-0.35 0.3];ytickBeta = [-0.35:0.1:0.3];
%10dots
ylimBeta= [2 4];ytickBeta = [2:0.2:4]; %comment out for baselined
%ylimBeta= [-1 0.3];ytickBeta = [-1:0.2:0.3];[] 

hf = figure(9)
hold on
for bin=1:1
    for c=1:2
        contra_beta_grp = squeeze(mean(mean(beta_side_all(:,c,[8  43],:,:),5),3));
       % ipsi_beta_grp = squeeze(mean(beta_side_all(:,c,25,:,:),5));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',ColsRT(c,:)}),
     %   h2((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
     %       {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[-100 1400],'ytick',ytickBeta,'xtick',[-100, 0:200:1400],'XTickLabel',xlabels_plot,'ylim',ylimBeta,'TickDir','out');%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Beta (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from stimulus onset (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
       % title('Contra Beta x Distractor xRT')
       line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
       line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
      % hss3((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
        %   'LineWidth',1.5,'LineStyle','--','Color',ColsRT(c,:));
       %hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
         %  'LineWidth',1.5,'LineStyle','--','Color',ColsRT(c,:));
    end
end
legend([h1.mainLine],dLabels,'FontSize',12,'Location','SouthEast')
legend boxoff
saveas(hf,['Figures\beta_stim_uB' namestring{onew} '.png']);
hf = figure(10)
hold on
for bin=1:1
    for c=1:2
        contra_beta_grp = squeeze(mean(mean(betar_side_all(:,c,[8 43],:,:),5),3));
       % ipsi_beta_grp = squeeze(mean(betar_side_all(:,c,25,:,:),5));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',ColsRT(c,:)}),
        %h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
        %    {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'linewidth',1.5,'Fontsmoothing','on','FontSize',16,'xlim',[-600 100],'ytick',ytickBeta,'xtick',-600:200:100,'ylim',ylimBeta,'TickDir','out');%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Beta (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time from Response (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
       % title('Contra Beta x Distractor xRT')
        line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','-','Color',ColsRT(c,:));
    end
end
legend([h1.mainLine],dLabels,'FontSize',12,'Location','SouthWest')
legend boxoff
saveas(hf,['Figures\beta_resp_uB' namestring{onew} '.png']);