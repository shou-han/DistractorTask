function [fast_slow] =  gen_bins(all_data, conds, no_of_bins)
% bins are associated with the trials
allD = all_data(conds);
[alldata_sort, indx] = sort(allD);
group_size = floor(length(alldata_sort)/no_of_bins);

fast_slow = zeros([1 no_of_bins]);
for bin = 1:no_of_bins
    indx((bin-1)*group_size+1);
    fast_slow(bin) = allD(indx((bin-1)*group_size+1));
end
end