function [onsetOut] = obtainOnset_final_subj( data,...
    search_window, consecutive_windows,...
    t, resplocked)
% extract the relevant data from search window
prestim_temp = find(t<search_window(1)); % so we can add it on after getting the data
dU = data(find(t>search_window(1) & t<search_window(2)),:); % time x trial
ts_temp=t(find(t>search_window(1) & t<search_window(2)));
% for CPP, best use resp-locked, for N2s use stimlocked
[mTime,mTrial] = size(dU);
if ~resplocked
    clear tstats ps allp05 onsetp05 
    for tt = 1:mTime
        clear ts sigDiffDatashow_
        sigDiffData = dU(tt,:);
        % do t-test to zero for the group of data
        [~,P,~,STATS] = ttest(sigDiffData);
        tstats(tt) = STATS.tstat;
        ps(tt) = P;
    end
    % find the last time x number of stats are significantly
    % positive
    allp05= find(ps<0.05 & tstats>0);
%     ps
%     tstats
%     %[ts_temp/1000;ps;tstats]
%     allp05

    onsetp05=[];
    for i = 1:length(allp05)
        if  (i+consecutive_windows-1)<=length(allp05)
            if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least 15 consecutive windows where p<.05
                onsetp05=allp05(i+consecutive_windows-2);
                break
            end
        else
            % no onset
            onsetp05=nan;
            break
        end
    end
    clear CPP_onset_ind
    length(prestim_temp);
    onsetp05;
    % get timepoint of min index.
    if ~isnan(onsetp05)
        CPP_onset_ind = onsetp05 + length(prestim_temp); % see above, this needs to be added to get the overall time with respect to t.
        onsetOut = t(CPP_onset_ind);
    else % onsetp05 is empty, no significant CPP.
        onsetOut = nan;
    end
end
% resplocked (for CPP)
if resplocked
    clear tstats ps allp05 onsetp05
    star = find(t==0);
    for tt = 1:star
        clear ts sigDiffData
        ts = star-tt+1;
        % look at allRT
        sigDiffData = dU(ts,:);
        % do t-test to zero for the group of data
        [~,P,~,STATS] = ttest(sigDiffData);
        tstats(tt) = STATS.tstat;
        ps(tt) = P;
    end
    % find the first time x number of stats are significantly
    % zero
    allp05= find(ps>0.05 & tstats>0);
    onsetp05=[];
    for i = 1:length(allp05)
        if  (i+consecutive_windows-1)<=length(allp05)
            if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least 15 consecutive windows where p>.05
                onsetp05=allp05(i+consecutive_windows-2);
                break
            end
        else %not 15 consecutive windows, can be lower, take half t
            onsetp05=nan;
            break
        end
    end
    clear CPP_onset_ind
    % get timepoint of min index.
    if ~isnan(onsetp05)
        CPP_onset_ind =  star-onsetp05+1;
        onsetOut= t(CPP_onset_ind); %note that this is resp-locked t
    else % onsetp05 is empty, no significant CPP.
        %         disp([allsubj,': bugger']) %AD48C has no CPP onset
        CPP_onset_ind = 1;
        onsetOut =  t(CPP_onset_ind);
    end
end
% plotting
% hold on
%
%plot(tr(CPP_onset_ind), data_used{1}(CPP_onset_ind,trial),'r+','Markersize',20)
%plot(t, data_used{1}(:,trial))
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

