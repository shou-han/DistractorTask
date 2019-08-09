function [dataOut] = fileCombiner( data,s,dataInPrev)

a = dataInPrev;
a.default(s,:,:,:,:) = data.default;
a.fast_RT(s,:,:,:,:)  = data.fast_RT;
a.slow_RT(s,:,:,:,:)  = data.slow_RT;
a.fast_e(s,:,:,:,:) = data.fast_e;
a.slow_e(s,:,:,:,:)  = data.slow_e;
a.correct(s,:,:,:,:)  = data.correct;
a.incorrect(s,:,:,:,:)  = data.incorrect;
a.high_rew(s,:,:,:,:)  = data.high_rew;
a.low_rew(s,:,:,:,:)  = data.low_rew;

dataOut= a';
end

