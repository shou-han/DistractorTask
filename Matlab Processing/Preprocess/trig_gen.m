function [trigAll ] = trig_gen( EEG, targcodes,max_response_time,fs)
%TRIG_GEN Summary of this function goes here
%   Detailed explanation goes here
numev = length(EEG.event)-1;

for i=1:numev
    trigAll.trigs(i)=EEG.event(i).type;
    trigAll.stimes(i)=round(EEG.event(i).latency);
end
targtrigs = [];
for n=1:length(trigAll.trigs)
    if any(targcodes(:)==trigAll.trigs(n))
        targtrigs = [targtrigs n];
    end
end

% Get rid of last trig, it was not recorded in the scores
trigAll.motion_on = targtrigs(1:end-1); 

% if trigAll.trigs(targtrigs(end))==trialCond(1)
% else
%     motion_on = targtrigs;
% end

trigAll.targtrigs = targtrigs;
trigAll.max_response_time = max_response_time;
end

