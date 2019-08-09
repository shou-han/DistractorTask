function measure = measure_gen(erp, allRT, ch_CPP)

        %% find CPP level at RT
        allRT_t = floor((allRT +700)/2+1);
        allRT_t(isnan(allRT_t)) = floor((1+700)/2+1);
        for chs=1:64
            erp_clk_temp(chs,:) = diag(squeeze(erp(chs,allRT_t,:))); %(s, chan, STFT_timers, targetside)
        end
        
        %% find CPP onset        
        %% Code adapted from Ger's Current Biology cpp code to pull out CPP onset latency:
        % Define CPP onset search window, from 0 to 1000ms
        CPP_search_t  = [10,1000];
        % Same window in samples
        CPP_search_ts  = [find(t==CPP_search_t(1)),find(t==CPP_search_t(2))];
        % Size of sliding window. This is in fact 1/4 of the search window in ms.
        % So 25 is 100ms. (25 samples x 2ms either side of a particular sample).
        max_search_window = 25;
        consecutive_windows=10;%10/15 works well for everybody else
        %% find onset CPP
        CPP_temp = squeeze(mean(erp(ch_CPP,:,[conds.default{s,c,:,side}]),1)); % time x trial

        % constrain the search window according to parameters above.
        CPPs(:,side) = squeeze(mean(CPP_temp(:,:),2)); % average across trial for plot later on, not used to find onsets.

        [CPP_side_onsets(s,c,side) , CPP_winm{side}]= obtainOnset(CPP_temp,CPP_search_t,t, ...
            max_search_window, 0.05, consecutive_windows, allsubj{s});
        
        
        % find CPP slope
        
        % find N2 onset
        
        % find N2 amplitude

end
plotting = 0;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_temp = '..\..\Data\';
fs = 500;
BL_erp = [-100,0];

% resp-locked erps
trs = [-.700*fs:fs*.100];

targcodes(1,:) = [101,103];
targcodes(2,:) = [102,104];
targcodes(3,:) = [105,107];
targcodes(4,:) = [106,108];
targcodes(5,:) = [109,111];
targcodes(6,:) = [110,112];

rtlim=[0.300 1.500];

ch_CPP = [31];
%% Use Current Source Density transformed erp? 1=yes, 0=no
CSD=0;
%% Start loop

for s= 1:18
    for c = 1:3%1:3
        load([path_temp subject_folder{s} '\' allsubj{s} '_' num2str(c) '_big_dots_erp'],'erp_LPF_8Hz','erp_LPF_35Hz','erp_LPF_35Hz_CSD','allRT','allrespLR','allTrig','allblock_count',...
            'BL_resp_artrej','ET_BL_resp_artrej')
        
        if CSD
            erp=double(erp_LPF_35Hz_CSD);
        else
            erp=double(erp_LPF_8Hz);
        end
        
        % Baseline erp
        baseline_erp = mean(erp(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
        erp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp
        
        disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(find(allTrig)))])
        
        erpr = zeros(size(erp,1),length(tr),size(erp,3));
        
        % check for valid reaction time in the required time
        validrlock = zeros(1,length(allRT)); % length of RTs.
        % for each trial, shift the erp signature to make RT = 0;
        for n=1:length(allRT)
            % shift time by the reaction time
            [blah,RTsamp] = min(abs(t*fs/1000-allRT(n)));
            if RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t) & allRT(n)>0 % is the RT larger than 1st stim RT point, smaller than 1.8seconds.
                erpr(:,:,n) = erp(:,RTsamp+trs,n);
                validrlock(n)=1;
            end
        end
        
        %% sorting out the data
        for side = 1:2
            for iti = 1:6
                % calcs the indices of the triggers for each appropriate trial type.
                conds.default{s,c,iti,side} = find(allTrig==targcodes(iti,side) & ...
                    allRT>rtlim(1)*fs & allRT<rtlim(2)*fs & BL_resp_artrej==1 & ET_BL_resp_artrej  & validrlock );
            end
        end
        
        for side = 1:2
        end
        
        %% CPP Onset Figures
        if plotting
            colors = {'b' 'r' 'g' 'm' 'c'};
            
            figure
            for side = 1:2
                plot(t,squeeze(CPPs(:,side)),'Color',colors{side},'LineWidth',2), hold on
                plot(10:2:1000,CPP_winm{side}','Color',colors{side},'LineWidth',2), hold on
                line([mean(CPP_side_onsets(s,c,side),1),mean(CPP_side_onsets(s,c,side),1)],ylim,'Color',colors{side},'LineWidth',1.5);
                line([0,0],ylim,'Color','k','LineWidth',1);
                line(xlim,[0,0],'Color','k','LineWidth',1);
            end
            pause
        end
        
    end
end

%% Extract Response Locked CPP slope"
%CPP build-up defined as the slope of a straight line fitted to the
%response-locked waveform at during "slope_timeframe_index" defined for
%each participant here:
clear slope_timeframe_index
for s=1:size(allsubj,2)
    for c=1:3
        slope_timeframe_index(s,c,2)=find(squeeze(mean(CPPr_side.default(s,c,:,:,:),5))==max(squeeze(mean(CPPr_side.default(s,c,:,find(tr<0),:),5))));%max amplitude index
        slope_timeframe_index(:,c,1)=slope_timeframe_index(:,c,2)-50;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    end
end
%Now find and save CPPr slope  for each subject
for s=1:size(allsubj,2)
    for c= 1:3
        for side = 1:2
            coef = polyfit(tr(slope_timeframe_index(s,c,1):slope_timeframe_index(s,c,2)),...
                squeeze(CPPr_side.default(s,c,:,slope_timeframe_index(s,c,1):slope_timeframe_index(s,c,2),side))'...
                ,1);% coef gives 2 coefficients fitting r = slope * x + intercept
            CPPr_slope(s,c,side)=coef(1);
            CPPr_int(s,c,side)=coef(2);
        end
    end
end

if plotting
    %%%Plot each individual participant's CPPr_slope with time-window varying
    %%%per participant
    for s=1:size(allsubj,2)
        pause(1)
        close all
        clear h
        for c= 1:3
            figure
            for side = 1:2
                h(side) = plot(tr,squeeze(CPPr_side.default(s,c,:,:,side)),'LineWidth',3,'LineStyle','-');hold on
                r = CPPr_slope(s,c,side) .* tr(slope_timeframe_index(s,c,1):slope_timeframe_index(s,c,2)) + CPPr_int(s,c,side); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
                plot(tr(slope_timeframe_index(s,c,1):slope_timeframe_index(s,c,2)), r,'Linewidth',2, 'LineStyle', ':');
            end
            
            set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
            ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('Time (ms)','FontName','Arial','FontSize',16)
            title([subject_folder{s}, ' CPP (resp-locked) by Hemifield'])
            line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
            line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
            legend(h,side_tags,'FontSize',16,'Location','NorthWest');
            pause(1)
        end
    end
end
%% Extract N2c and N2i peak latency :
for s = 1:size(allsubj,2)
    for c = 1:3
        for side=1:2
            %         to use each participant's average N2c/i to get their peak latency index:
            N2c=squeeze(N2c_side.default(s,c,:, :, side));
            N2i=squeeze(N2i_side.default(s,c,:, :, side));
            N2c_peak_amp_index_t(s,c,side)=t(N2c==min(N2c(find(t==150):find(t==400))));%Find max peak latency for N2c in ms
            N2i_peak_amp_index_t(s,c,side)=t(N2i==min(N2i(find(t==200):find(t==450))));%Find max peak latency for N2i in ms
        end
    end
end
%N2 Latency:
N2cN2i_latency_ByTargetSide = [N2c_peak_amp_index_t,N2i_peak_amp_index_t]; %(LeftTargetN2c_latency, RightTargetN2c_latency, LeftTargetN2i_latency, RightTargetN2i_latency)
if plotting
    %%%Plot each individual participant's CPPr_slope with time-window varying
    %%%per participant
    for s=1:size(allsubj,2)
        pause(1)
        close all
        clear h
        for c= 1:3
            figure
            for side = 1:2
                N2c=squeeze(N2c_side.default(s,c,:, :, side));
                h(side) = plot(t,squeeze(N2c_side.default(s,c,:,:,side)),'LineWidth',3,'LineStyle','-');hold on
                r = min(N2c(find(t==150):find(t==400)));
                hold on
                line([N2c_peak_amp_index_t(s,c,side)-window, N2c_peak_amp_index_t(s,c,side)+window],...
                    [r r],'Color','k','Linewidth',2, 'LineStyle', '-');
            end
            
            set(gca,'FontSize',16,'xlim',[-100,1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
            ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('Time (ms)','FontName','Arial','FontSize',16)
            title([subject_folder{s}, ' N2c'])
            line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
            line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
            legend(h,side_tags,'FontSize',16,'Location','NorthWest');
            pause
        end
    end
end
%% Extract N2c and N2i Amplitude :
window=25; %this is the time (in samples) each side of the peak latency (so it's 50ms each side of peak latency - so a 100ms window)
for c= 1:3
    clear N2c N2i
    N2c = squeeze(mean(mean(N2c_side.default(:,c,:,:,:),1),5)); % time
    N2i = squeeze(mean(mean(N2i_side.default(:,c,:,:,:),1),5)); % time
    N2c_peak_amp_index(c)=find(N2c==min(N2c(find(t==150):find(t==450))));%Find Left target max peak latency for N2c
    N2i_peak_amp_index(c)=find(N2i==min(N2i(find(t==150):find(t==450))));%Find Left target max peak latency for N2i
end

for s = 1:size(allsubj,2)
    for c = 1:3
        for side=1:2
            max_peak_N2c(s,c,side)=squeeze(min(N2c_side.default(s,c,:,N2c_peak_amp_index(c)-window:N2c_peak_amp_index(c)+window, side),[],4));
            max_peak_N2i(s,c,side)=squeeze(min(N2i_side.default(s,c,:,N2i_peak_amp_index(c)-window:N2i_peak_amp_index(c)+window, side),[],4));
        end
    end
end

N2cN2i_amp_ByTargetSide_ParticipantLevel = [max_peak_N2c,max_peak_N2i]; %(LeftTargetN2c, RightTargetN2c, LeftTargetN2i, RightTargetN2i)
if plotting
    %%%Plot each individual participant's CPPr_slope with time-window varying
    %%%per participant
    for s=1:size(allsubj,2)
        pause(1)
        close all
        clear h
        for c= 1:3
            figure
            for side = 1:2
                h(side) = plot(t,squeeze(N2c_side.default(s,c,:,:,side)),'LineWidth',3,'LineStyle','-');hold on
                r = max_peak_N2c(s,c,side);
                hold on
                line([t(N2c_peak_amp_index(c)-window), t(N2c_peak_amp_index(c)+window)],...
                    [r r],'Color','k','Linewidth',2, 'LineStyle', '-');
            end
            
            set(gca,'FontSize',16,'xlim',[-100,1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
            ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('Time (ms)','FontName','Arial','FontSize',16)
            title([subject_folder{s}, ' N2c'])
            line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
            line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
            legend(h,side_tags,'FontSize',16,'Location','NorthWest');
            pause(1)
        end
    end
    
    %%
    %%%Plot each individual participant's N2i with time-window varying
    %%%per participant
    for s=1:size(allsubj,2)
        pause(1)
        close all
        clear h
        for c= 1:3
            figure
            for side = 1:2
                h(side) = plot(t,squeeze(N2i_side.default(s,c,:,:,side)),'LineWidth',3,'LineStyle','-');hold on
                r = max_peak_N2i(s,c,side);
                hold on
                line([t(N2i_peak_amp_index(c)-window), t(N2i_peak_amp_index(c)+window)],...
                    [r r],'Color','k','Linewidth',2, 'LineStyle', '-');
            end
            
            set(gca,'FontSize',16,'xlim',[-100,1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
            ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('Time (ms)','FontName','Arial','FontSize',16)
            title([subject_folder{s}, ' N2i'])
            line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
            line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
            legend(h,side_tags,'FontSize',16,'Location','NorthWest');
            pause(1)
        end
    end
end
%% Make participant level matrix for export into SPSS or R
for c = 1:3
    participant_level(:,1:2)=squeeze(max_peak_N2c(:,c,:)); %N2c amplitude (LeftTarget, RightTarget)
    participant_level(:,3:4)=squeeze(max_peak_N2i(:,c,:)); %N2i amplitude (LeftTarget, RightTarget)
    participant_level(:,5:6)=squeeze(N2c_peak_amp_index_t(:,c,:)); %N2c latency (LeftTarget, RightTarget)
    participant_level(:,7:8)=squeeze(N2i_peak_amp_index_t(:,c,:)); %N2i latency (LeftTarget, RightTarget)
    participant_level(:,9:10)=squeeze(CPP_side_onsets(:,c,:)); %CPP onset (LeftTarget, RightTarget)
    participant_level(:,11:12)=squeeze(CPPr_slope(:,c,:)); %response locked CPP slope (LeftTarget, RightTarget)
    % open participant_level

    csvwrite (['participant_level_matrix' num2str(c) '.csv'],participant_level)

    subject_folder=subject_folder';
    cell2csv (['IDs' num2str(c) '.csv'],subject_folder)
end