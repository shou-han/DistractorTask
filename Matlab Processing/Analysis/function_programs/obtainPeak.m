function [peakMag, peakTime] = obtainPeak( data, ...
    search_window, consecutive_windows,...
    t, resplocked)
% extract the relevant data from search window
prestim_temp = find(t<search_window(1)); % so we can add it on after getting the data
dU = data(find(t>search_window(1) & t<search_window(2)),:); % time x trial
dt = t(find(t>search_window(1) & t<search_window(2)));

% for CPP, best use resp-locked, for N2s use stimlocked

[mTime,~] = size(dU);
if ~resplocked
    clear tstats ps allp05 onsetp05 ts timeData
    for tt = 1:mTime
        % find the highest time
        timeData(tt) = mean(dU(tt,:));
    end
    peakData =max(timeData);
    pT_temp = find(timeData == peakData,1);
    % keep going until significantly different
    for tt = 1:pT_temp
        % do t-test to zero for the group of data
        [~,P,~,STATS] = ttest2(dU(pT_temp-tt+1,:), dU(pT_temp,:));
        tstats(tt) = STATS.tstat;
        ps(tt) = P;
    end
    % find the first time x number of stats are significantly
    % different
    allp05= find(ps<0.05);
   % ps
    onsetp05=[];
    for i = 1:max(length(allp05),1)
        if  (i+consecutive_windows-1)<=length(allp05)
            if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least n consecutive windows where p<.05
                onsetp05=allp05(i+consecutive_windows-2);
                break
            end
        else
            % all the same
            onsetp05=pT_temp;
            break
        end
    end
   % pT_temp
   % onsetp05
    onsetp05 = pT_temp - onsetp05+1; %convert to actual time
   % onsetp05
   % t(pT_temp+length(prestim_temp))
   % t(onsetp05+length(prestim_temp))
    clear CPP_onset_ind
    % get timepoint of min index.
        CPP_onset_ind = onsetp05 + length(prestim_temp); % see above, this needs to be added to get the overall time with respect to t.
        onsetOut = t(CPP_onset_ind);
end
% resplocked (for CPP)

if resplocked
    clear tstats ps allp05 onsetp05 ts timeData
    % find mean 
    for tt = 1:mTime
        % find the highest time
        timeData(tt) = mean(dU(tt,:));
    end
    peakData =max(timeData);
    pT_temp = find(timeData == peakData,1);
    for tt = 1:pT_temp
        clear ts sigDiffData
        ts = pT_temp-tt+1;
        % look at allRT
        sigDiffData = peakData - dU(ts,:);
        % do t-test to zero for the group of data
        [~,P,~,STATS] = ttest(sigDiffData);
        tstats(tt) = STATS.tstat;
        ps(tt) = P;
    end
    % find the first time x number of stats are significantly
    % zero
    allp05= find(ps<0.05);
    onsetp05=[];
    for i = 1:max(length(allp05),1)
        if  (i+consecutive_windows-1)<=length(allp05)
            if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least 15 consecutive windows where p>.05
                onsetp05=allp05(i+consecutive_windows-2);
                break
            end
        else %not 15 consecutive windows, can be lower, take half t
            onsetp05=pT_temp;
            break
        end
    end
    clear CPP_onset_ind
    % get timepoint of min index.
    CPP_onset_ind =  pT_temp-onsetp05+1;
    onsetOut= t(CPP_onset_ind); %note that this is resp-locked t
end
% plotting
% hold on
% %
% figure
% hold on
% 
% %plot(t(CPP_onset_ind), mean(data(CPP_onset_ind,:)),'r+','Markersize',20)
% plot(t, mean(data,2))
% plot(dt, mean(dU,2),'r')
% plot(dt,timeData,'y')
% plot(t(CPP_onset_ind), mean(data((CPP_onset_ind),:),2),'r+','Markersize',20)
% plot(t(pT_temp+length(prestim_temp)), mean(data((pT_temp+length(prestim_temp)),:),2),'r+','Markersize',20)
% pause

peakMag = peakData; 
peakTime = t(pT_temp+length(prestim_temp));
%plot(tr(RT_used(trial)), data_used{1}(RT_used(trial),trial),'r+','Markersize',20);

%CPP_search_t,t,maxT,max_search_window, ptS, , allsubj, ch_CPP,c,side, totalOnsetOut, totalAmpIndx, totalAmp, resplocked)

% data_used{2} = squeeze(mean(data(ch_CPP,:,[conds.fast_RT{c,:,side}]),1)); % time x trial
% data_used{3} = squeeze(mean(data(ch_CPP,:,[conds.slow_RT{c,:,side}]),1)); % time x trial
% data_used{4} = squeeze(mean(data(ch_CPP,:,[conds.fast_e{c,:,side}]),1)); % time x trial
% data_used{5} = squeeze(mean(data(ch_CPP,:,[conds.slow_e{c,:,side}]),1)); % time x trial
% data_used{6} = squeeze(mean(data(ch_CPP,:,[conds.accurate{c,:,side}]),1)); % time x trial
% data_used{7} = squeeze(mean(data(ch_CPP,:,[conds.inaccurate{c,:,side}]),1)); % time x trial
% data_used{8} = squeeze(mean(data(ch_CPP,:,[conds.high_rew{c,:,side}]),1)); % time x trial
% data_used{9} = squeeze(mean(data(ch_CPP,:,[conds.low_rew{c,:,side}]),1)); % time x trial
% totalOnsetOut.fast_RT(s,c,side) = onsetOut(2);
% totalOnsetOut.slow_RT(s,c,side) = onsetOut(3);
% totalOnsetOut.fast_e(s,c,side) = onsetOut(4);
% totalOnsetOut.slow_e(s,c,side) = onsetOut(5);
% totalOnsetOut.correct(s,c,side) = onsetOut(6);
% totalOnsetOut.incorrect(s,c,side) = onsetOut(7);
% totalOnsetOut.high_rew(s,c,side) = onsetOut(8);
% totalOnsetOut.low_rew(s,c,side) = onsetOut(9);

%totalAmpIndx.default(c,:,side) = ampIndex{1}.trials;
% totalAmpIndx.fast_RT{s,c,side} = ampIndex{2};
% totalAmpIndx.slow_RT{s,c,side} = ampIndex{3};
% totalAmpIndx.fast_e{s,c,side} = ampIndex{4};
% totalAmpIndx.slow_e{s,c,side} = ampIndex{5};
% totalAmpIndx.correct{s,c,side} = ampIndex{6};
% totalAmpIndx.incorrect{s,c,side} = ampIndex{7};
% totalAmpIndx.high_rew{s,c,side} = ampIndex{8};
% totalAmpIndx.low_rew{s,c,side} = ampIndex{9};

%totalAmp.default(c,:,side) = maxAmp{1}.trials;
% totalAmp.fast_RT{s,c,side} = maxAmp{2};
% totalAmp.slow_RT{s,c,side} = maxAmp{3};
% totalAmp.fast_e{s,c,side} = maxAmp{4};
% totalAmp.slow_e{s,c,side} = maxAmp{5};
% totalAmp.correct{s,c,side} = maxAmp{6};
% totalAmp.incorrect{s,c,side} = maxAmp{7};
% totalAmp.high_rew{s,c,side} = maxAmp{8};
% totalAmp.low_rew{s,c,side} = maxAmp{9};
end

