function [energyLevel] = obtainEnergy( data_used, t, RT_used,resplocked,flag)
%obtainbetaclick( data_used, t,RT_used, resplocked,'sum'/'mean')

% find beta at click
if ~isequal( size(data_used,2), size(RT_used,2))
    error('erp and reaction times trial numbers are not equal')
end
%% baseline it to zero
%baseline = mean(data_used{1}(351-increment_window_size:351 + increment_window_size,:));
%data_used{1} = data_used{1}-repmat(baseline, [size(data_used{1},1),1]);

% whether resplocked or not, will always search for the onset as the point
% closest to the time of the click which is not significantly different
% from zero.


% not resplocked 
%let RT_used be 0 and then the previous times as positive time
% keep on going until either find the non-zero or until the time gets to
% RT_used
if ~resplocked
    for trial = 1:size(data_used,2)
        t_zero_index= find(t<0,1,'last');
        [RT_closest, ts]= min(abs(t-RT_used(trial)));
        %ts = RT_index+t_zero_index;
        % look at allRT
        sigDiffData = data_used(t_zero_index:ts,trial);
        energyLevelTemp(trial) = norm(sigDiffData);
    end
    
end

% resplocked
if resplocked
    t_zero_index= find(t<0,1,'last');
    for trial = 1:size(data_used,2)
            ts = t_zero_index;
            % look at allRT
            %size(data_used{1})
            sigDiffData = data_used(1:ts,trial);
            energyLevelTemp(trial) = norm(sigDiffData);
    end
end

if strcmp(flag,'sum')
energyLevel = sum(energyLevelTemp);
end
if strcmp(flag,'mean')
energyLevel = mean(energyLevelTemp);
end
% plotting
% hold on
% trial = 34;
% plot(tr(451), CPP_click(flds).trials(trial),'r+','Markersize',20)
% plot(tr, data_used{1}(:,trial))
% plot(tr(451),data_used{1}(451,trial),'b+','Markersize',20)

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

