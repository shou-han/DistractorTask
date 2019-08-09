function [totalOnsetOut, totalAmpIndx, totalAmp] = obtainOnset_final( data, conds, CPP_search_t,t,maxT,max_search_window, ptS, consecutive_windows, allsubj, ch_CPP,s,c,side, totalOnsetOut, totalAmpIndx, totalAmp, resplocked)
%% find onset CPP
data_used{1} = squeeze(mean(data(ch_CPP,:,[conds.default{s,c,:,side}]),1)); % time x trial
data_used{2} = squeeze(mean(data(ch_CPP,:,[conds.fast_RT{s,c,:,side}]),1)); % time x trial
data_used{3} = squeeze(mean(data(ch_CPP,:,[conds.slow_RT{s,c,:,side}]),1)); % time x trial
data_used{4} = squeeze(mean(data(ch_CPP,:,[conds.fast_e{s,c,:,side}]),1)); % time x trial
data_used{5} = squeeze(mean(data(ch_CPP,:,[conds.slow_e{s,c,:,side}]),1)); % time x trial
data_used{6} = squeeze(mean(data(ch_CPP,:,[conds.accurate{s,c,:,side}]),1)); % time x trial
data_used{7} = squeeze(mean(data(ch_CPP,:,[conds.inaccurate{s,c,:,side}]),1)); % time x trial
data_used{8} = squeeze(mean(data(ch_CPP,:,[conds.high_rew{s,c,:,side}]),1)); % time x trial
data_used{9} = squeeze(mean(data(ch_CPP,:,[conds.low_rew{s,c,:,side}]),1)); % time x trial

tstart = CPP_search_t(1);
if resplocked
    tstart= t(1);
end
for flds = 1:9
    clearvars -except flds data_used  conds CPP_search_t t max_search_window maxT ptS consecutive_windows allsubj ch_CPP s c side onsetOut totalOnsetOut ampIndex totalAmpIndx tstart totalAmp maxAmp 
    % find(t>=CPP_search_t(1) & t<=CPP_search_t(2))
    % constrain the search window according to parameters above.
    
    dU = data_used{flds};
    
    if size(dU,1) == 1
        dU = dU';
    end
    
    if isempty(dU)
        dU = zeros*data_used{1};
    end
    CPP_temp = squeeze(dU(find(t>=CPP_search_t(1) & t<=CPP_search_t(2)),:));
    prestim_temp = find(t<CPP_search_t(1)); % so we can add it on after getting max peak.
    % only looking at the times from 0 to 1000ms
    % moving average the
    win_mean_inds = 1:1:size(CPP_temp,1);
    for trial = 1:size(CPP_temp,2)
        win_mean(:,trial) = movmean(CPP_temp(win_mean_inds,trial),max_search_window);
    end
    
    % do t-test to zero across the smoothed trials.
    for tt = 1:size(win_mean,1)
        [~,P,~,STATS] = ttest(win_mean(tt,:));
        tstats(tt) = STATS.tstat;
        ps(tt) = P;
    end
    
    % find when the 10 tstats are significantly positive for all
    % trials
    clear allp05
    allp05= find(ps<ptS & tstats>0);
    onsetp05=[];
    for i = 1:length(allp05)
        if  (i+consecutive_windows-1)<=length(allp05)
            if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least 10 consecutive windows where p<.05
                onsetp05=allp05(i);
                break
            end
        else
            onsetp05=allp05(i);
            break
        end
    end
    
    % get timepoint of min index.
    if ~isempty(onsetp05)
        onset_ind = win_mean_inds(onsetp05); %from 0 to 1000 window
        CPP_onset_ind = onset_ind + length(prestim_temp);
        onsetOut(flds) = t(CPP_onset_ind);
    else % onsetp05 is empty, no significant CPP.
        
        %         disp([allsubj,': bugger']) %AD48C has no CPP onset
        CPP_onset_ind = length(prestim_temp);
        onsetOut(flds) =  t(CPP_onset_ind);
    end

    dtemp = data_used{flds};
    if size(dtemp,1) == 1
        dtemp = dtemp';
    end

    if ~isempty(dtemp)
%         onsetOut(flds)
%         max(t)
temp2 = mean(dtemp(find(t>onsetOut(flds),1):find(t==maxT),:),2);
            maxAmp{flds} = max(temp2);
            ampIndex{flds} = find( temp2 == maxAmp{flds},1)+CPP_onset_ind;
            clear temp2
    else % no trials

        ampIndex{flds} = CPP_onset_ind; %no amplitude
        maxAmp{flds} = 0;
%              disp([allsubj,': bugger'])
%              pause
    end
    clear dtemp
    
    
%     % find the maximum amplitude of CPP
%        close(figure(1101))
%     figure(1101)
%     hold on
%     ampM = floor(mean([ampIndex{flds}])); maxA = (mean([maxAmp{flds}]));
%     plot(t,mean(dU,2))
%     line([onsetOut(flds) onsetOut(flds)],[-0.1 0.1],'LineWidth',1.5,'LineStyle','--');
%     line([t(ampM) t(ampM)],[-0.1 0.1],'LineWidth',1.5,'LineStyle','--');
%         line([t(ampM-2), t(ampM+2)],[maxA maxA],'LineWidth',1.5,'LineStyle','--');
%     pause
end
if isempty(onsetOut(1) )
     disp([allsubj,': bugger'])
end
totalOnsetOut.default(s,c,side) = onsetOut(1);
totalOnsetOut.fast_RT(s,c,side) = onsetOut(2);
totalOnsetOut.slow_RT(s,c,side) = onsetOut(3);
totalOnsetOut.fast_e(s,c,side) = onsetOut(4);
totalOnsetOut.slow_e(s,c,side) = onsetOut(5);
totalOnsetOut.correct(s,c,side) = onsetOut(6);
totalOnsetOut.incorrect(s,c,side) = onsetOut(7);
totalOnsetOut.high_rew(s,c,side) = onsetOut(8);
totalOnsetOut.low_rew(s,c,side) = onsetOut(9);

totalAmpIndx.default{s,c,side} = ampIndex{1};
totalAmpIndx.fast_RT{s,c,side} = ampIndex{2};
totalAmpIndx.slow_RT{s,c,side} = ampIndex{3};
totalAmpIndx.fast_e{s,c,side} = ampIndex{4};
totalAmpIndx.slow_e{s,c,side} = ampIndex{5};
totalAmpIndx.correct{s,c,side} = ampIndex{6};
totalAmpIndx.incorrect{s,c,side} = ampIndex{7};
totalAmpIndx.high_rew{s,c,side} = ampIndex{8};
totalAmpIndx.low_rew{s,c,side} = ampIndex{9};

totalAmp.default{s,c,side} = maxAmp{1};
totalAmp.fast_RT{s,c,side} = maxAmp{2};
totalAmp.slow_RT{s,c,side} = maxAmp{3};
totalAmp.fast_e{s,c,side} = maxAmp{4};
totalAmp.slow_e{s,c,side} = maxAmp{5};
totalAmp.correct{s,c,side} = maxAmp{6};
totalAmp.incorrect{s,c,side} = maxAmp{7};
totalAmp.high_rew{s,c,side} = maxAmp{8};
totalAmp.low_rew{s,c,side} = maxAmp{9};
end

