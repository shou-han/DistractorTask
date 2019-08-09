function [CPPmean] = obtainCPP_click_final_subj( data,  allRT, increment_window_size, t,resplocked)
%% find CPP at time of click
data_used=data;
%% baseline it to zero
%baseline = mean(data_used{1}(351-increment_window_size:351 + increment_window_size,:));
%data_used{1} = data_used{1}-repmat(baseline, [size(data_used{1},1),1]);
RT_used = allRT; % sample time at x trial - used only if not resp_locked
%
% whether resplocked or not, will always search for the onset as the point
% closest to the time of the click which is not significantly different
% from zero.


% not resplocked 
%let RT_used be 0 and then the previous times as positive time
% keep on going until either find the non-zero or until the time gets to
% RT_used
if ~resplocked
    for trial = 1:size(data_used,2)
            ts = find(t>RT_used(trial),1);
            % look at allRT
            sigDiffData = data_used(ts-increment_window_size:ts + increment_window_size,trial);
            CPP_click(trial) = mean(sigDiffData);
    end
    
end

% resplocked
if resplocked
    for trial = 1:size(data_used,2)
            ts = find(t==0);
            % look at allRT
            %size(data_used{1})
            sigDiffData = data_used(ts-increment_window_size:ts + increment_window_size,trial);
            CPP_click(trial) = mean(sigDiffData);
    end
end
% plotting
% hold on
% trial = 34;
% plot(tr(451), CPP_click(flds).trials(trial),'r+','Markersize',20)
% plot(tr, data_used{1}(:,trial))
% plot(tr(451),data_used{1}(451,trial),'b+','Markersize',20)
CPPmean = mean(CPP_click);
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

