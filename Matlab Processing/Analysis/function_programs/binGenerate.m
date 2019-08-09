function youBin = binGenerate(ynoBin, mediator, conds, no_of_bins)
           %separate into conditions
            RT_temp = []; conds_bin_temp = [];
            RT_temp = [allRTdiff([conddiff{s,cs,:,side}])];
            %bin it
            [RT_temp,indx] = sort(RT_temp);
            cond_bin_temp = [conddiff{s,cs,:,side}];
            conds_bin_temp  = cond_bin_temp(indx);
            group_size = floor(length(RT_temp)/no_bins);
            rem_size = rem(length(RT_temp),no_bins);
            
            %choose the smallest bin index
            binch = find(bins_count_diff==min(bins_count_diff),1);
            %choose the smallest bin index
            %             if rem_size>0
            %                 binch = randint(1,1,[1 no_of_bins]);
            %             else
            %                 binch = 1;
            %             end
            
            for bin =  [1:binch-1]
                condsdiffBin{s,cs,side,bin} = conds_bin_temp((bin-1)*group_size+1:bin*group_size);
                bins_count_diff(bin) = bins_count_diff(bin)+group_size;
            end
            
            condsdiffBin{s,cs,side,binch} = conds_bin_temp((binch-1)*group_size+1:(binch)*group_size+rem_size); %account for the remainders
            bins_count_diff(binch) = bins_count_diff(binch)+group_size+rem_size;
            
            for bin = [binch+1:no_bins]
                condsdiffBin{s,cs,side,bin} = conds_bin_temp((bin-1)*group_size+rem_size+1:bin*group_size+rem_size);
                bins_count_diff(bin) = bins_count_diff(bin)+group_size;
            end
end