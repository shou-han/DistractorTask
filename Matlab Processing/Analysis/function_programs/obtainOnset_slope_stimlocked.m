function [totalOnsetOut1, totalOnsetOut2] = obtainOnset_slope( data,search_window, CPPr_onset, tr)
% slope is defined as the gradient at the  half way point between the onset and the  final
dU= data;
onsetTime = min(search_window(1), CPPr_onset);
%% baseline it to zero
%baseline = mean(data_used{1}(351-increment_window_size:351 + increment_window_size,:));
%data_used{1} = data_used{1}-repmat(baseline, [size(data_used{1},1),1]);
%
% slope now tries to find the maximum velocity given the position, so sort
% of derivative, but using all the trial distributions instead of two
% points. 
clear cpp_slope_data_tot time_data_tot 

st = find(tr> onsetTime,1);
ed = find(tr==50,1);
cpp_slope_data_tot = dU(st:ed,:)';
tr_tot = tr(st:ed);
[m,n] = size(cpp_slope_data_tot);
% if n<70
%     tstep = floor(n/2); %take half the points as gradient
% else
%     tstep =50; %20  points
% end
CPPr_slope = [];
clear cpp_slope_data time_data
cpp_slope_data = reshape(cpp_slope_data_tot',n*m,1)';
time_data = repmat(tr_tot,[1 m]);
%size(cpp_slope_data)
%size(time_data)
coef = fitlm(time_data,cpp_slope_data,'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
CPPr_slope=coef.Coefficients.Estimate(2);
CPPr_int = coef.Coefficients.Estimate(1);



totalOnsetOut1 = CPPr_slope;
totalOnsetOut2 = CPPr_int;

end

