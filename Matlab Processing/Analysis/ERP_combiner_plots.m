clear all
close all
set(0,'DefaultFigureVisible','on')
t=1;
addpath(genpath('/scratch/yn70/ShouHan/Distractors/CSDtoolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/eeglab13_6_5b/'));
addpath('function_programs/');
%eeglab
%ssj = [1:20];
ssj=[1:18 20:21];
namestring ={'150dots','10dots'};
%ssj=[5];
%ssj = [1:3 7:13 16:22 24];
%ssj=[1:13 15:16 18:22 ];%5 14
%ssj = [8 10 12:20 21 22:26 29 ];
addpath('function_programs/');
onew=2;
oldnew={'', 'old'};
%ssj = [1:5 8:26 28:29];
chn = 13;
CSD = 1;
if CSD
    yCPP = [-10,35]; ytickCPP =[-10:5:35];
    yN2pc = [-30,20]; ytickN2pc = [-30:5:20];
    yN2c= [-30,20]; ytickN2c=[-30:5:20];
    yN2i=[-30,20];ytickN2i=[-30:5:20];
    ylim = [-30 30]; ylimBeta=[-2 0.5];ytickBeta=[-2:0.3:0.5];
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
%combine them
for s = 1:size(ssj,2)
    sT = ssj(s);
    if CSD
        load(['Data/ERPs/group_plots_erp_diff_CSD_1binsStats_8Hz_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','t','tr','allBins','no_of_bins',...
            'beta_side','beta_r_side','beta_contra','betar_contra','STFT_time','STFT_timer','RT_side',...
            'betaslope','alpha_side','alpha_asym_side','alpha_t','alpha_tr')
    else
        load(['Data/ERPs/group_plots_erp_diff_1bin_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'RTs'...
            ,'chanlocs','t','tr','allBins','no_of_bins',...
            'beta_side','beta_r_side','beta_contra','betar_contra','STFT_time','STFT_timer','RT_side',...
            'betaslope','alpha_side','alpha_asym_side','alpha_t','alpha_tr')
    end
    
    allBins = allBins;
    %allBinStats;
    % conds=cat(1,conds,condsdiff);
    for bin=1:no_of_bins
        for d=1:length(cs)
            for cong=1
                c = cs(d,:);
                %all_stuff.RTs{s,c} = ([RTs{1,c,:,:}]);
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
red         = rgb('DeepPink');
blue        = [0    0.7461    1.0000];
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
%%
% ERP
for bin=1
for c=1:2
    dataplottopor(:,:,c) = squeeze(mean(mean(mean(ERPr_sides(:,c,:,:,:),1),2),5));
end
end
figure(500)
plottopo(dataplottopor,'chanlocs',chanlocs,'limits',[tr(1) tr(size(tr,2)) ...
    min(min(min(dataplottopor))) max(max(max(dataplottopor)))], ...
    'title',['ERPr hand and distractors'],'showleg','on','ydir',1)
%%
%CPPR
for c=1:2
    dataplottoporbeta(:,:,c) = squeeze(mean(mean(mean(betar_side_all(:,c,:,:,:),1),2),5));
end
figure(902)
plottopo(dataplottoporbeta,'chanlocs',chanlocs,'limits',[STFT_timer(1) STFT_timer(end) ...
    min(min(min(dataplottoporbeta))) max(max(max(dataplottoporbeta)))], ...
    'title',['ERPr hand and distractors'],'showleg','on','ydir',1)
%topoplot
%% Topology Plots
%% ERP scalp topo
if CSD
    maplim = [-15 15];
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

for c=2
    %N2pc
    hf1 = figure(1)
    subplot(1,1,1)
    plot_mean(left_hemi) = squeeze(mean(mean(mean(mean(ERP_sides(:,c,left_hemi,find(t>=280 & t<350),:)-...
        ERP_sides(:,1,right_hemi,find(t>=280 & t<350),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',left_hemi);
    %N2c/N2i
    hf2=figure(2)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,c,:,find(t>=350 & t<450),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',1:64);
    %N2c/N2i
    hf3=figure(3)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,c,:,find(t>=450 & t<550),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        maplim, ...
        'electrodes','off','plotchans',1:64);
    %CPP
    hf4=figure(4)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,c,:,find(t>=550 & t<650),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
     [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
    %BetaWeird
    hf5=figure (5)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(beta_side_all(:,c,:,find(STFT_time>=600 & STFT_time<800),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean),max(plot_mean)], ...
        'electrodes','off','plotchans',1:64);
    %Beta
    hf6=figure (6)
    subplot(1,1,1)
    plot_mean = squeeze(mean(mean(mean(mean(betar_side_all(:,c,:,find(STFT_timer>=-20 & STFT_timer<20),:),1),2),4),5));
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(-0.9),max(0.7)], ...
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
return
%% Topoplots
for c=1:2
    dataplottopor(:,:,c) = squeeze(mean(mean(mean(ERPr_sides(:,c,:,:,:),1),2),5));
end
figure(500)
plottopo(dataplottopor,'chanlocs',chanlocs,'limits',[tr(1) tr(end) ...
    -50 50], ...
    'title',['ERPr hand and distractors'],'showleg','on','ydir',1)
%% Beta Topoplots
clear dataplottopo
for c=1:2
    dataplottopo(:,:,c) = squeeze(mean(mean(mean(beta_side_all(:,c,:,:,:),1),2),5));
end
figure(500)
plottopo(dataplottopo,'chanlocs',chanlocs,'limits',[STFT_time(1) STFT_time(end) ...
    -0.5 0.5], ...
    'title',['beta hand and distractors'],'showleg','on','ydir',1)


%% plottopo
bhf = figure (903)
plot_mean = squeeze(mean(mean(mean(mean(betar_side_all(:,2,:,find(STFT_timer>=-50 & STFT_timer<50),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [-2 1], ...
    'electrodes','labels','plotchans',1:64);
colorbar

bhf = figure (904)
plot_mean = squeeze(mean(mean(mean(mean(betar_side_all(:,2,:,find(STFT_timer>=-50 & STFT_timer<50),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [-2 1], ...
    'electrodes','labels','plotchans',1:64);
colorbar

%% RT Behaviour
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
      set(gca,'FontSize',16,'xlim',[400,1600],'xtick',[400:200:1600],'ylim',[0,500],'ytick',[0:150:450]);%,'ylim',[-1.5,0.5]);
    h = histc(de_temp, [min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)]);
    hs(cs) =plot([min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)],h,'Color',ColsRT(c,:),'Linewidth',1.5);
    hb(cs) = line([mean(meanRTs{c}) mean(meanRTs{c})],[0 450],'Color',ColsRT(c,:),'Linewidth',1.5,'LineStyle','--');
    
end
ylabel({'Frequency across';' all Participants'},'FontSize',16,'fontweight','bold')
xlabel('Response Time (ms)','FontSize',16,'fontweight','bold')
legend(hs,{'DA','DP'},'FontSize',12,'Location','NorthEast','fontweight','bold')
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
ylabel({'Response Time (ms)'},'FontSize',16,'fontweight','bold')
axis([0 3 300 1200]);
xticks([1 2]);xticklabels({'DA','DP'});
set(gca,'Fontsize',20)
yticks(300:300:1200);
legend(hs,{'Distractor Absent','Distractor Present'},'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\RT_bar_graphs' namestring{onew} '.png']);
%% N2 graphs
%shadedErrorBar
%N2pc
clear hf hss1
hf = figure(1)
for bin=1:1
    for c=1:2
        N2pc_c = squeeze(mean(mean(N2pc_sides(:,c,:,:),4),2));
        meanN2pc = mean(N2pc_c,1);
        stdN2pc = std(N2pc_c,1)/sqrt(size(N2pc_c,1));
        hold on
        (c-1)*2+bin
        hss1((bin-1)*2+c) = shadedErrorBar(t,meanN2pc,stdN2pc,'lineprops',{'Color',[dpurple alphaRT((bin-1)*2+c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[],'ylim',yN2pc,'ytick',ytickN2pc);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            ylim,'Color',[dpurple alphaRT((bin-1)*2+c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2pc Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        %xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16','fontweight','bold')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine, hss1(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2pcxdist' namestring{onew} '.png']);
% N2c 
clear hss3
hf = figure(3)
for bin = 1:1
    hold on
    for c=1:2
        N2c_c_grp = squeeze(mean(mean(ERP_sides(:,c,45,:,:,:),2),5));
        meanN2c = mean(N2c_c_grp,1);
        stdN2c = std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[dred alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[],...
            'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2c,'Color',[dred alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
       % xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title(' N2c x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2cxdist' namestring{onew} '.png']);
% N2i 
clear hss3
hf = figure(2)
for bin = 1:1
    hold on
    for c=1:2
        N2i_c_grp = squeeze(mean(mean(ERP_sides(:,c,51,:,:,:),2),5));
        meanN2i = mean(N2i_c_grp,1);
        stdN2i = std(N2i_c_grp,1)/sqrt(size(N2i_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[dblue alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',yN2i,'ytick',ytickN2i);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2i,'Color',[dblue alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2i Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
       % title(' N2c x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2ixdist' namestring{onew} '.png']);


%% CPP_CPPr 
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
        hss1((bin-1)*2+c) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1400],'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        %hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
         %   ylim,'Color',[dgreen alphaRT((bin-1)*2+c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('CPP Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine hss1(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
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
        hss11((bin-1)*2+c) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('from response (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
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
        hss11((bin-1)*2+c) = shadedErrorBar(alpha_t,meanalpha,stdalpha,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0 16],'ytick',[0:2:16]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('from stimulus (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
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
        hss11((bin-1)*2+c) = shadedErrorBar(alpha_t,meanalpha,stdalpha,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0.1 0.1],'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-0.5 0.5],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('from stimulus (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\alpha_sides_asym' namestring{onew} '.png']);
%%
clear hss11
hf = figure(6)
for bin = 1:1
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPPslope_sides(:,c,:,:),2),4));
        CPP_cs = squeeze(mean(CPP_c));
        CPPstd = squeeze(std(CPP_c))/sqrt(size(CPP_c,1));
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(t,CPP_cs,CPPstd,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',-100:200:1400,...
            'ylim',[-0.5 0.5],'ytick',[-0.5:0.05:0.5]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-1 1],'Color',[Cols(c+3*(bin-1),:)],'LineWidth',1.5,'LineStyle','--');
        ylabel('CPP Slope (\muVolts/sec)','FontName','Arial','FontSize',16)
       % xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
       % title('CPP slopes')
        line([0 0],ylim,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPslope.png']);

clear hss11
hf = figure(7)
for bin = 1:1
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPPrslope_sides(:,c,:,:),4),2));
        CPP_cs = squeeze(mean(CPP_c));
        CPPstd = squeeze(std(CPP_c))/sqrt(size(CPP_c,1));
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(tr,CPP_cs,CPPstd,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'})
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',-600:100:100,...
            'ylim',[-0.5 0.5],'ytick',[]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        %ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        %xlabel('from response (ms)','FontName','Arial','FontSize',16)
        %title('CPPr Slope')
        line([0 0],ylim,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
%legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPrSlope.png']);
%% beta
clear h1 h2 hss4 hf hss3
% contra beta and ipsi beta x distractors resp
hf = figure(9)
hold on
for bin=1:1
    for c=1:2
        contra_beta_grp = squeeze(mean(beta_side_all(:,c,8,:,:),5));
       % ipsi_beta_grp = squeeze(mean(beta_side_all(:,c,25,:,:),5));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
     %   h2((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
     %       {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'FontSize',16,'xlim',[-100 1400],'ytick',[-2:0.3:0.5],'xtick',[-100, 0:200:1400],'ylim',[-2 0.5]);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
      %  xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title('Contra Beta x Distractor xRT')
       line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
       line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
       hss3((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
           'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]);
       hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
           'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]);
    end
end
legend([h1.mainLine],{'DA','DP'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\beta_stim' namestring{onew} '.png']);
hf = figure(10)
hold on
for bin=1:1
    for c=1:2
        contra_beta_grp = squeeze(mean(betar_side_all(:,c,8,:,:),5));
       % ipsi_beta_grp = squeeze(mean(betar_side_all(:,c,25,:,:),5));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
        %h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
        %    {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'FontSize',16,'xlim',[-600 100],'ytick',[-2:0.3:0.5],'xtick',-600:200:100,'ylim',[-2 0.5]);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title('Contra Beta x Distractor xRT')
        line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','-','Color',[RTCols alphaRT(c)]);
    end
end
legend([h1.mainLine],{'DA','DP'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\beta_resp' namestring{onew} '.png']);

clear h1 h2
hf = figure(11)
for bin = 1:1
    hold on
    for c=1:2
        beta_c = squeeze(mean(mean(beta_c_slope_all(:,c,:,:),2),4));
      %  beta_i = squeeze(mean(mean(beta_i_slope_all(:,c,:,:),2),4));        
        beta_cs = squeeze(mean(beta_c));%beta_is = squeeze(mean(beta_i));
        beta_cstd = squeeze(std(beta_c));%beta_istd = squeeze(mean(beta_i));
        hold on
        h1((bin-1)*2+c) = shadedErrorBar(STFT_time,beta_cs,beta_cstd,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
      %  h2((bin-1)*2+c) = shadedErrorBar(STFT_time,beta_is,beta_istd,'lineprops',...
       %     {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}), 
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[],...
            'ylim',[-0.1 0.1],'ytick',[-0.1:0.1:0.1]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-1 1],'Color',[ColsBeta(1+2*(bin-1),:)],'LineWidth',1.5,'LineStyle','--');
        ylabel('Beta Slope (\muVolts/sec)','FontName','Arial','FontSize',16)
       % xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
       % title('CPP slopes')
        line([0 0],[-0.1 0.1],'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([h1(1).mainLine,h1(2).mainLine],{'DA','DP'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\betaslope' namestring{onew} '.png']);
clear h1 h2
hf = figure(12)
for bin = 1:1
    hold on
    for c=1:2
        beta_c = squeeze(mean(mean(betar_c_slope_all(:,c,:,:),2),4));
      %  beta_i = squeeze(mean(mean(beta_i_slope_all(:,c,:,:),2),4));        
        beta_cs = squeeze(mean(beta_c));%beta_is = squeeze(mean(beta_i));
        beta_cstd = squeeze(std(beta_c));%beta_istd = squeeze(mean(beta_i));
        hold on
        h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,beta_cs,beta_cstd,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
      %  h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,beta_is,beta_istd,'lineprops',...
       %     {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}), 
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[],...
            'ylim',[-0.1 0.1],'ytick',[-0.1:0.1:0.1]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-1 1],'Color',[ColsBeta(1+2*(bin-1),:)],'LineWidth',1.5,'LineStyle','--');
        ylabel('Beta Slope (\muVolts/sec)','FontName','Arial','FontSize',16)
       % xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
       % title('CPP slopes')
        line([0 0],[-0.1 0.1],'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([h1(1).mainLine,h1(2).mainLine],{'DA','DP'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\betarslope' namestring{onew} '.png']);
%% Topology Plots
%% ERP scalp topo
if CSD
    maplim = [-15 15];
    mapRlim = [-20 20];
    mapDlim = [-10 10];  
else
    maplim = [-1.5 1.5];
    mapRlim = [-2 2];
    mapDlim = [-1 1];
end
%%%%%%%%%%%%%%%%%DA
bhf = figure (902)
left_hemi = [1 33 34 4 3 37 36 5 38 6 39 7 9 41 8 40 10 42 11 43 12 15 45 14 44 46 47 16];
right_hemi = [32 63 62 31 30 60 61 27 59 28 58 29 26 56 25 57 21 55 22 54 23 20 51 19 52 50 49 18];
centre_chans = [35 48 2 13 17 64 24 53];

plot_mean = zeros([64,1]);
plot_mean(left_hemi) = squeeze(mean(mean(mean(mean(ERP_sides(:,1,left_hemi,find(t>=280 & t<420),:)-...
    ERP_sides(:,1,right_hemi,find(t>=280 & t<420),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',left_hemi);
saveas(bhf,['Figures\N2pctopo_DA.png'])

bhf = figure (903)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,1,:,find(t>=280 & t<320),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
     [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\N2ctopo_DA.png'])

bhf = figure (904)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,1,:,find(t>=380 & t<420),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    maplim, ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\N2itopo_DA.png'])

bhf = figure (905)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,1,:,find(t>=700 & t<800),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
     [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\CPPtopo_DA.png'])

bhf = figure (906)
plot_mean = squeeze(mean(mean(mean(mean(beta_side_all(:,1,:,find(STFT_time>=700 & STFT_time<900),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\Betatopo_DA.png'])

bhf = figure (907)
plot_mean = squeeze(mean(mean(mean(mean(betar_side_all(:,1,:,find(STFT_timer>=-100 & STFT_timer<100),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\Betatopo_DA.png'])
%%%%%%%%%%%%%%%%%DP
bhf = figure (908)
plot_mean = zeros([64,1]);
plot_mean(left_hemi) = squeeze(mean(mean(mean(mean(ERP_sides(:,2,left_hemi,find(t>=280 & t<420),:)-...
    ERP_sides(:,1,right_hemi,find(t>=280 & t<420),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',left_hemi);
saveas(bhf,['Figures\N2pctopo_DP.png'])

bhf = figure (909)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,2,:,find(t>=280 & t<320),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\N2ctopo_DP.png'])

bhf = figure (910)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,2,:,find(t>=380 & t<420),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    maplim, ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\N2itopo_DP.png'])

bhf = figure (911)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,2,:,find(t>=680 & t<720),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\CPPtopo_DP.png'])

bhf = figure (912)
plot_mean = squeeze(mean(mean(mean(mean(beta_side_all(:,2,:,find(STFT_time>=700 & STFT_time<900),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\Betatopo_DP.png'])

bhf = figure (913)
plot_mean = squeeze(mean(mean(mean(mean(betar_side_all(:,2,:,find(STFT_timer>=-50 & STFT_timer<50),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [-0.13,max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\Betatopo_DP.png'])
%%%%%% 
%% RT Bins

%% Average Response Times
clear Rt_sj_times RT_group_times
for bin=1:no_of_bins
    for sj = 1:size(allstuff_a.RTs,1)
        
        for c = 1:2
            Rt_sj_times(sj,c,bin) = mean([allstuff_a.RTs{sj,c,bin}]);
        end
    end
    RT_group_times(:,:,bin) = Rt_sj_times(:,:,bin)';
end
%% CPP x RT Bins
hf = figure(41)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPP_sides(:,c,:,bin),4),2));
        meanCPP = mean(CPP_c,1);
        stdCPP = std(CPP_c,1)/sqrt(size(CPP_c,1));             
        hold on
        hss1((bin-1)*1+c) =  shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',[Cols(1+(bin-1),:) alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100 0:200:1400],...
            'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*1+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yCPP,'Color',[Cols(1+(bin-1),:) alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('CPP Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss1(1).mainLine, hss1(2).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPxRT.png']);
% CPPr x RT Bins
clear hss11 hss4
hf = figure(51)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        CPPr_c = squeeze(mean(mean(CPPr_sides(:,c,:,bin),4),2));
        meanCPPr = mean(CPPr_c,1);
        stdCPPr = std(CPPr_c,1)/sqrt(size(CPPr_c,1));
        hold on
        hss11((bin-1)*1+c) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',[Cols(1+(bin-1),:) alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:0],...
            'ylim',yCPP,'ytick',ytickCPP);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yCPP,'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('from response (ms)','FontName','Arial','FontSize',16)
        %title('CPPr x RT Bins')
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPrxRT.png']);
%%
clear hss11
hf = figure(61)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPPslope_sides(:,c,:,bin),2),4));
        CPP_cs = squeeze(mean(CPP_c));
        CPPstd = squeeze(std(CPP_c))/sqrt(size(CPP_c,1));
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(t,CPP_cs,CPPstd,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',[-0.1 0.1],'ytick',[-0.1:0.05:0.1]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-1 1],'Color',[Cols(c+3*(bin-1),:)],'LineWidth',1.5,'LineStyle','--');
        ylabel('CPP Slope (\muVolts/sec)','FontName','Arial','FontSize',16)
        xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
       % title('CPP slopes')
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPslopexrt.png']);

clear hss11
hf = figure(71)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        CPP_c = squeeze(mean(mean(CPPrslope_sides(:,c,:,bin),4),2));
        CPP_cs = squeeze(mean(CPP_c));
        CPPstd = squeeze(std(CPP_c))/sqrt(size(CPP_c,1));
        hold on
        hss11((bin-1)*2+c) = shadedErrorBar(tr,CPP_cs,CPPstd,'lineprops',{'Color',[Cols(c+3*(bin-1),:)],'LineWidth',3,'LineStyle','-'})
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[-600:200:100],...
            'ylim',[-0.1 0.1],'ytick',[]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-10 30],'Color',[RTCols alphaRT(c)],'LineWidth',1.5,'LineStyle','-');
        %ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('from response (ms)','FontName','Arial','FontSize',16)
        %title('CPPr Slope')
        line([0 0],yCPP,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
%legend([hss11(1).mainLine, hss11(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\CPPrSlopexrt.png']);

%% N2xRT Bins
%N2pc
clear hss3 hss4
hf = figure(4)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        N2pc_c_grp = squeeze(mean(mean(N2pc_sides(:,c,:,bin),2),6));
        meanN2pc = mean(N2pc_c_grp,1);
        stdN2pc = std(N2pc_c_grp,1)/sqrt(size(N2pc_c_grp,1));        
        hold on
        hss3((bin-1)*1+c) = shadedErrorBar(t,meanN2pc,stdN2pc,'lineprops',{'Color',[Cols(7+(bin-1),:) alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[],...
            'ylim',yN2pc,'ytick',[]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*1+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2pc,'Color',[Cols(7+(bin-1),:) alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        %ylabel('N2pc Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        %xlabel('Time (ms)','FontName','Arial','FontSize',16)
        %title(' N2pc x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2pcxRT.png']);
% N2cxRT
hf = figure(5)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        N2c_c_grp = squeeze(mean(mean(ERP_sides(:,c,45,:,bin,:),2),6));
        meanN2c = mean(N2c_c_grp,1);
        stdN2c = std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));        
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*1+c) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[Cols(2+3*(bin-1),:) alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[],...
            'ylim',yN2c,'ytick',[]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*1+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2c,'Color',[Cols(2+3*(bin-1),:) alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        %ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        %xlabel('Time (ms)','FontName','Arial','FontSize',16)
        %title(' N2i x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2cxRT.png']);
% N2ixRT
hf = figure(6)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        N2i_c_grp = squeeze(mean(mean(ERP_sides(:,c,51,:,bin,:),2),6));
        meanN2i = mean(N2i_c_grp,1);
        stdN2i = std(N2i_c_grp,1)/sqrt(size(N2i_c_grp,1));          
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*1+c) = shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[Cols(3+3*(bin-1),:) alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',yN2i,'ytick',[]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*1+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2i,'Color',[Cols(3+3*(bin-1),:) alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        %ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
        %title(' N2i x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\N2ixRT.png']);

%% Beta xRT
clear hss3 hss4
hf = figure(11)
hold on
for bin=1:no_of_bins
    for c=1:2
        contra_beta_grp = squeeze(mean(beta_contra_all(:,c,:,bin),4));
        %ipsi_beta_grp = squeeze(mean(beta_ipsi_all(:,c,:,bin),4));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
        %h2((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
        %   {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'FontSize',16,'xlim',[-100 1400],'ytick',ytickBeta,'xtick',[],'ylim',ylimBeta);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Contra Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title('Contra Beta x Distractor xRT')
       line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
       line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
       hss3((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
           'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]);
     %  hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
      %     'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]);
    end
end
legend([h1(2).mainLine,h1(4).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','SouthWest')
saveas(hf,'Figures\beta_stim_RT.png');
hf = figure(12)
hold on
for bin=1:no_of_bins
    for c=1:2
        contra_beta_grp = squeeze(mean(betar_contra_all(:,c,:,bin),4));
       % ipsi_beta_grp = squeeze(mean(betar_ipsi_all(:,c,:,bin),4));
        mean_contra_beta = mean(contra_beta_grp);%mean_ipsi_beta = mean(ipsi_beta_grp);
        std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);%std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
        h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_contra_beta,std_contra_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
%        h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
 %           {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'FontSize',16,'xlim',[-600 100],'ytick',[],'xtick',[],'ylim',ylimBeta);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
    %    ylabel('Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
      %  xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title('Contra Beta x Distractor xRT')
        line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','-','Color',[RTCols alphaRT(c)]);
    end
end
legend([h1(2).mainLine,h1(4).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','SouthWest')
saveas(hf,'Figures\beta_resp_RT.png');
%
hf = figure(13)
hold on
for bin=1:no_of_bins
    for c=1:2
   %     contra_beta_grp = squeeze(mean(beta_contra_all(:,c,:,bin),4));
        ipsi_beta_grp = squeeze(mean(beta_ipsi_all(:,c,:,bin),4));
     %   mean_contra_beta = mean(contra_beta_grp);
        mean_ipsi_beta = mean(ipsi_beta_grp);
     %   std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);
        std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
   %     h1((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_contra_beta,std_contra_beta,'lineprops',...
    %        {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
        h2((bin-1)*2+c) = shadedErrorBar(STFT_time,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'FontSize',16,'xlim',[-100 1400],'ytick',ytickBeta,'xtick',[-100, 0:200:1400],'ylim',ylimBeta);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
        ylabel('Ipsi Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
      %  xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title('Contra Beta x Distractor xRT')
       line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
       line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    %   hss3((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
     %      'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]);
       hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
           'LineWidth',1.5,'LineStyle','--','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]);
    end
end
legend([h2(2).mainLine,h2(4).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','SouthWest')
saveas(hf,'Figures\beta_ipsi_stim_RT.png');
%ipsi
hf = figure(14)
hold on
for bin=1:no_of_bins
    for c=1:2
       % contra_beta_grp = squeeze(mean(betar_contra_all(:,c,:,bin),4));
        ipsi_beta_grp = squeeze(mean(betar_ipsi_all(:,c,:,bin),4));
        %mean_contra_beta = mean(contra_beta_grp);
        mean_ipsi_beta = mean(ipsi_beta_grp);
       % std_contra_beta = std(contra_beta_grp)/size(contra_beta_grp,1);
       std_ipsi_beta = std(ipsi_beta_grp)/size(ipsi_beta_grp,1);
  %      h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_contra_beta,std_contra_beta,'lineprops',...
   %         {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
        h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,mean_ipsi_beta,std_ipsi_beta,'lineprops',...
           {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}),      
        set(gca,'FontSize',16,'xlim',[-600 100],'ytick',[],'xtick',[-600:200:100],'ylim',ylimBeta);%'ylim',[0.5 0.8]);%,'ylim',[0.5,0.8]);
    %    ylabel('Beta Amplitude (\muVolts)','FontName','Arial','FontSize',16)
      %  xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title('Contra Beta x Distractor xRT')
        line([0,0],ylimBeta,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        hss4((bin-1)*2+c)=line([mean(RT_group_times(c,:,bin),2), mean(RT_group_times(c,:,bin),2)],ylimBeta,...
            'LineWidth',1.5,'LineStyle','-','Color',[RTCols alphaRT(c)]);
    end
end
legend([h2(2).mainLine,h2(4).mainLine],{'Fast RT','Slow RT'},'FontSize',12,'Location','SouthWest')
saveas(hf,'Figures\beta_ipsi_resp_RT.png');

clear h1 h2
hf = figure(15)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        beta_c = squeeze(mean(mean(beta_c_slope_all(:,c,:,bin),2),4));
      %  beta_i = squeeze(mean(mean(beta_i_slope_all(:,c,:,bin),2),4));        
        beta_cs = squeeze(mean(beta_c));%beta_is = squeeze(mean(beta_i));
        beta_cstd = squeeze(0*std(beta_c));%beta_istd = squeeze(mean(beta_i));
        hold on
        h1((bin-1)*2+c) = shadedErrorBar(STFT_time,beta_cs,beta_cstd,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
      %  h2((bin-1)*2+c) = shadedErrorBar(STFT_time,beta_is,beta_istd,'lineprops',...
       %     {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}), 
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[],...
            'ylim',[-0.5 0.5],'ytick',[-0.5:0.1:0.5]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-1 1],'Color',[ColsBeta(1+2*(bin-1),:)],'LineWidth',1.5,'LineStyle','--');
        ylabel('Beta Slope (\muVolts/sec)','FontName','Arial','FontSize',16)
       % xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
       % title('CPP slopes')
        line([0 0],ylim,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([h1(1).mainLine,h1(2).mainLine],{'DA','DP'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\betaslope_RT.png']);
clear h1 h2
hf = figure(16)
for bin = 1:no_of_bins
    hold on
    for c=1:2
        beta_c = squeeze(mean(mean(betar_c_slope_all(:,c,:,bin),2),4));
      %  beta_i = squeeze(mean(mean(beta_i_slope_all(:,c,:,bin),2),4));        
        beta_cs = squeeze(mean(beta_c));%beta_is = squeeze(mean(beta_i));
        beta_cstd = squeeze(std(beta_c));%beta_istd = squeeze(mean(beta_i));
        hold on
        h1((bin-1)*2+c) = shadedErrorBar(STFT_timer,beta_cs,beta_cstd,'lineprops',...
            {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(1+(bin-1)*2,:) alphaRT(c)]}),
      %  h2((bin-1)*2+c) = shadedErrorBar(STFT_timer,beta_is,beta_istd,'lineprops',...
       %     {'LineWidth',3,'LineStyle','-','Color',[ColsBeta(2+(bin-1)*2,:) alphaRT(c)]}), 
        set(gca,'FontSize',16,'xlim',[-600,100],'xtick',[],...
            'ylim',[-0.5 0.5],'ytick',[-0.5:0.1:0.5]);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-1 1],'Color',[ColsBeta(1+2*(bin-1),:)],'LineWidth',1.5,'LineStyle','--');
        ylabel('Beta Slope (\muVolts/sec)','FontName','Arial','FontSize',16)
       % xlabel('Time from evidence (ms)','FontName','Arial','FontSize',16)
       % title('CPP slopes')
        line([0 0],ylim,'Color','k','LineWidth',3.0,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([h1(1).mainLine,h1(2).mainLine],{'DA','DP'},'FontSize',12,'Location','SouthWest')
saveas(hf,['Figures\betarslope_RT.png']);
%% presentation graphs
%% RT Behaviour
clear hs hb
    for c = 1:2
        temp  =[];
        for s = 1:length(ssj)
            temp = [temp [allstuff_a.RTs{s,c,:}]];
        end
        meanRTs{c} = temp;
        clear temp 
    end
hf = figure(131)
hold on
cs=0;
for c = [1 2]
    cs=cs+1;
    clear de_temp
    de_temp = [allstuff_a.RTs{:,c,:}];
      set(gca,'FontSize',16,'xlim',[400,1600],'xtick',[400:200:1600],'ylim',[0,500],'ytick',[0:150:450]);%,'ylim',[-1.5,0.5]);
    h = histc(de_temp, [min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)]);
    hs(cs) =plot([min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)],h,'Color',ColsRT(c,:),'Linewidth',1.5);
    hb(cs) = line([mean(meanRTs{c}) mean(meanRTs{c})],[0 450],'Color',ColsRT(c,:),'Linewidth',1.5,'LineStyle','--');
    
end
ylabel({'Frequency across';' all Participants'},'FontSize',16,'fontweight','bold')
xlabel('Response Time (ms)','FontSize',16,'fontweight','bold')
legend(hs,{'DA','DP'},'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\Presentation\RT_bins_present.png']);

left_hemi = [1 33 34 4 3 37 36 5 38 6 39 7 9 41 8 40 10 42 11 43 12 15 45 14 44 46 47 16];
right_hemi = [32 63 62 31 30 60 61 27 59 28 58 29 26 56 25 57 21 55 22 54 23 20 51 19 52 50 49 18];
centre_chans = [35 48 2 13 17 64 24 53];

plot_mean = zeros([64,1]);

bhf = figure (132)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,1,:,find(t>=280 & t<320),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
     [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\Presentation\N2ctopo_DA_present.png'])
bhf = figure (133)
plot_mean = squeeze(mean(mean(mean(mean(ERP_sides(:,2,:,find(t>=280 & t<320),:),1),2),4),5));
topoplot(plot_mean,chanlocs,'maplimits', ...
     [min(plot_mean),max(plot_mean)], ...
    'electrodes','off','plotchans',1:64);
saveas(bhf,['Figures\Presentation\N2ctopo_DP_present.png'])
% N2c 
clear hss3
hf = figure(134)
for bin = 1:1
    hold on
    for c=1:2
        N2c_c_grp = squeeze(mean(mean(ERP_sides(:,c,45,:,:,:),2),5));
        meanN2c = mean(N2c_c_grp,1);
        stdN2c = std(N2c_c_grp,1)/sqrt(size(N2c_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',[dred alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',yN2c,'ytick',ytickN2c);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            yN2c,'Color',[dred alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2c Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
       % title(' N2c x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\Presentation\N2cxdist.png']);
% N2i 
clear hss3
hf = figure(135)
for bin = 1:1
    hold on
    for c=1:2
        N2i_c_grp = squeeze(mean(mean(ERP_sides(:,c,51,:,:,:),2),5));
        meanN2i = mean(N2i_c_grp,1);
        stdN2i = std(N2i_c_grp,1)/sqrt(size(N2i_c_grp,1));
        hold on
        %hss1((bin-1)*2+c) = plot(t,CPP_c,'Color',[Cols(1,:) alphas((bin-1)*2+c)],'LineWidth',3,'LineStyle','-');
        hss3((bin-1)*2+c) =shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',[dblue alphaRT(c)],'LineWidth',3,'LineStyle','-'});
        set(gca,'FontSize',16,'xlim',[-100,1400],'xtick',[-100,0:200:1200],...
            'ylim',yN2i,'ytick',ytickN2i);%,'ylim',[-1.5,0.5]);
        hss4((bin-1)*2+c)=line([mean(mean(RT_group_times(c,:,bin),1),2) mean(mean(RT_group_times(c,:,bin),1),2)],...
            [-21 10],'Color',[dblue alphaRT(c)],'LineWidth',1.5,'LineStyle','--');
        ylabel('N2i Amplitude (\muVolts)','FontName','Arial','FontSize',16,'fontweight','bold')
        xlabel('Time (ms)','FontName','Arial','FontSize',16,'fontweight','bold')
       % title(' N2c x RT Bins')
        line([0 0],ylim,'Color','k','LineWidth',3,'LineStyle','-');
        line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');
        %     legend(h,side_tags, ...
        %         'FontSize',16,'Location','NorthWest');
    end
end
legend([hss3(1).mainLine, hss3(2).mainLine],{'DA','DP'},'FontSize',12,'Location','NorthWest')
saveas(hf,['Figures\Presentation\N2ixdist.png']);