function [totalOnsetOut] = obtainOnset_old_subj_tl( data, ch_CPP, conds, winsmth, consecutive_windows, totalOnsetOut,...
    c,side,t)
%% find onset CPP
CPP_search_t = [0 1000];
data_used{1} = squeeze(mean(data(ch_CPP,:,[conds.default{c,:,side}]),1)); % time x trial
% data_used{2} = squeeze(mean(data(ch_CPP,:,[conds.fast_RT{c,:,side}]),1)); % time x trial
% data_used{3} = squeeze(mean(data(ch_CPP,:,[conds.slow_RT{c,:,side}]),1)); % time x trial
% data_used{4} = squeeze(mean(data(ch_CPP,:,[conds.fast_e{c,:,side}]),1)); % time x trial
% data_used{5} = squeeze(mean(data(ch_CPP,:,[conds.slow_e{c,:,side}]),1)); % time x trial
% data_used{6} = squeeze(mean(data(ch_CPP,:,[conds.accurate{c,:,side}]),1)); % time x trial
% data_used{7} = squeeze(mean(data(ch_CPP,:,[conds.inaccurate{c,:,side}]),1)); % time x trial
% data_used{8} = squeeze(mean(data(ch_CPP,:,[conds.high_rew{c,:,side}]),1)); % time x trial
% data_used{9} = squeeze(mean(data(ch_CPP,:,[conds.low_rew{c,:,side}]),1)); % time x trial
%%
%plot(t,dU);
%hold on
%plot(t(floor(allRT(1))+351), dU(floor(allRT(1))+351),'r+','Markersize',15)
%%
for flds = 1:1
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
    prestim_temp = find(t<0); % so we can add it on after getting max peak.
    % only looking at the times from 0 to 1000ms
    % moving average the
    win_mean_inds = 1:1:size(CPP_temp,1);
    for trial = 1:size(CPP_temp,2)
        win_mean(:,trial) = movmean(CPP_temp(win_mean_inds,trial),winsmth);
        %win_mean(:,trial) = abs(win_mean(:,trial));
    end
        % do t-test to zero across the smoothed trials.
    for tt = 1:size(win_mean,1)
        [~,P,~,STATS] = ttest(win_mean(tt,:));
        tstats(tt) = STATS.tstat;
        ps(tt) = P;
    end
    
    % find when the 10 tstats are significantly zero for all
    % trials
    clear allp05
    allp05= find(ps<0.05 & tstats>0);
    onsetp05=[];
    for i = 1:length(allp05)
        if  (i+consecutive_windows-1)<=length(allp05)
            if allp05(i+consecutive_windows-1)-allp05(i)==consecutive_windows-1 %if there is at least 10 consecutive windows where p>.05
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
    
%    plot(t(CPP_onset_ind), mean(dU(CPP_onset_ind,:),2),'r+','Markersize',20)
%    hold on
%    plot(t, mean(dU(:,:),2))
%    pause
    
    totalOnsetOut.default(c,side) = onsetOut(1);
    % totalOnsetOut.fast_RT(s,c,side) = onsetOut(2);
    % totalOnsetOut.slow_RT(s,c,side) = onsetOut(3);
    % totalOnsetOut.fast_e(s,c,side) = onsetOut(4);
    % totalOnsetOut.slow_e(s,c,side) = onsetOut(5);
    % totalOnsetOut.correct(s,c,side) = onsetOut(6);
    % totalOnsetOut.incorrect(s,c,side) = onsetOut(7);
    % totalOnsetOut.high_rew(s,c,side) = onsetOut(8);
    % totalOnsetOut.low_rew(s,c,side) = onsetOut(9);

end

