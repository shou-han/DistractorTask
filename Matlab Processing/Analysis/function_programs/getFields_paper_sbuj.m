function out  = getFields_paper_sbuj(data,chn,conds,c,side, in, options)
% 1 = eeg
% 2 = behaviour
out = in;

if options ~=3 && options < 5
    
    kks{1} = [conds.default{c,:,side}];
    kks{2} = [conds.fast_RT{c,:,side}];
    kks{3} = [conds.slow_RT{c,:,side}];
    kks{4} = [conds.fast_e{c,:,side}];
    kks{5} = [conds.slow_e{c,:,side}];
    kks{6} = [conds.accurate{c,:,side}];
    kks{7} = [conds.inaccurate{c,:,side}];
    kks{8} = [conds.high_rew{c,:,side}];
    kks{9} = [conds.low_rew{c,:,side}];
    
    for iii = 1:9
        if options == 1
            if ~isempty(kks{iii})
                outdata{iii} = squeeze(mean(data(:,:,kks{iii}),3));
            else
                outdata{iii} = 0*squeeze(mean(data(:,:,1),3));
            end
        end
        if options == 2
            if ~isempty(kks{iii})
                outdata{iii} = data(kks{iii});
            else
                outdata{iii} = 0*data(1);
            end
        end
        if options == 4
            if ~isempty(kks{iii})
                outdata{iii} = squeeze(mean(data(:,kks{iii}),2));
            else
                outdata{iii} = 0*data(:,1);
            end
        end
    end
    
    if options == 2
        out.default{c,side} = outdata{1};
        out.fast_RT{c,side}= outdata{2};
        out.slow_RT{c,side} = outdata{3};
        out.fast_e{c,side} = outdata{4};
        out.slow_e{c,side} = outdata{5};
        out.correct{c,side} = outdata{6};
        out.incorrect{c,side} = outdata{7};
        out.high_rew{c,side} = outdata{8};
        out.low_rew{c,side} = outdata{9};
    end
    if options ==1
        out.default(c,:,:,side) = outdata{1};
        out.fast_RT(c,:,:,side)= outdata{2};
        out.slow_RT(c,:,:,side) = outdata{3};
        out.fast_e(c,:,:,side) = outdata{4};
        out.slow_e(c,:,:,side) = outdata{5};
        out.correct(c,:,:,side) = outdata{6};
        out.incorrect(c,:,:,side) = outdata{7};
        out.high_rew(c,:,:,side) = outdata{8};
        out.low_rew(c,:,:,side) = outdata{9};
    end
    if options ==4
        
        out.default(c,:,side) = outdata{1};
        out.fast_RT(c,:,:,side)= outdata{2};
        out.slow_RT(c,:,:,side) = outdata{3};
        out.fast_e(c,:,side) = outdata{4};
        out.slow_e(c,:,side) = outdata{5};
        out.correct(c,:,side) = outdata{6};
        out.incorrect(c,:,side) = outdata{7};
        out.high_rew(c,:,side) = outdata{8};
        out.low_rew(c,:,side) = outdata{9};
    end
end

if options ==3 || options == 7
    
    kks{1} = [conds.defaultE{c,:,side}];
    kks{2} = [conds.fast_e{c,:,side}];
    kks{3} = [conds.slow_e{c,:,side}];
    kks{4} = [conds.accurateE{c,:,side}];
    kks{5} = [conds.inaccurateE{c,:,side}];
    kks{6} = [conds.high_rewE{c,:,side}];
    kks{7} = [conds.low_rewE{c,:,side}];
    
    for iii = 1:7
        if options == 3
            if ~isempty(kks{iii})
                outdata{iii} = squeeze(mean(data(:,:,kks{iii}),3));
            else
                outdata{iii} = 0*squeeze(mean(data(:,:,1),3));
            end
        end
        if options == 7
            if ~isempty(kks{iii})
                outdata{iii} = data(kks{iii});
            else
                outdata{iii} = 0*data(1);
            end
        end
        
    end
    
    if options == 3
    out.default(c,:,:,side) = outdata{1};
    out.fast_e(c,:,:,side) = outdata{2};
    out.slow_e(c,:,:,side) = outdata{3};
    out.correct(c,:,:,side) = outdata{4};
    out.incorrect(c,:,:,side) = outdata{5};
    out.high_rew(c,:,:,side) = outdata{6};
    out.low_rew(c,:,:,side) = outdata{7};
    end
    
    if options == 7
        out.default{c,side} = outdata{1};
        out.fast_e{c,side} = outdata{2};
        out.slow_e{c,side} = outdata{3};
        out.correct{c,side} = outdata{4};
        out.incorrect{c,side} = outdata{5};
        out.high_rew{c,side} = outdata{6};
        out.low_rew{c,side} = outdata{7};
    end
    
    
end

%% get all the data as one long graph
if options ==5
    
    kks{1} = [conds.default{c,:,side}];
    kks{2} = [conds.fast_RT{c,:,side}];
    kks{3} = [conds.slow_RT{c,:,side}];
    kks{4} = [conds.fast_e{c,:,side}];
    kks{5} = [conds.slow_e{c,:,side}];
    kks{6} = [conds.accurate{c,:,side}];
    kks{7} = [conds.inaccurate{c,:,side}];
    kks{8} = [conds.high_rew{c,:,side}];
    kks{9} = [conds.low_rew{c,:,side}];
    
    indata{1} = out.default{c,side};
    indata{2} = out.fast_RT{c,side};
    indata{3} = out.slow_RT{c,side};
    indata{4} = out.fast_e{c,side};
    indata{5} = out.slow_e{c,side};
    indata{6} = out.correct{c,side};
    indata{7} = out.incorrect{c,side};
    indata{8} = out.high_rew{c,side};
    indata{9} = out.low_rew{c,side};
    
    for iii = 1:9
        if ~isempty(kks{iii})
            if size(kks{iii},2)>1
                outdata{iii} = [indata{iii} squeeze(data(chn,:,kks{iii}))];
            else
                outdata{iii} = [indata{iii}];
            end
        else
            outdata{iii} = [indata{iii}];
        end
    end
    
    out.default{c,side} = outdata{1};
    out.fast_RT{c,side}= outdata{2};
    out.slow_RT{c,side}= outdata{3};
    out.fast_e{c,side} = outdata{4};
    out.slow_e{c,side}= outdata{5};
    out.correct{c,side}= outdata{6};
    out.incorrect{c,side}= outdata{7};
    out.high_rew{c,side} = outdata{8};
    out.low_rew{c,side} = outdata{9};
end

if options ==6
    
    kks{1} = [conds.defaultE{c,:,side}];
    kks{2} = [conds.fast_e{c,:,side}];
    kks{3} = [conds.slow_e{c,:,side}];
    kks{4} = [conds.accurateE{c,:,side}];
    kks{5} = [conds.inaccurateE{c,:,side}];
    kks{6} = [conds.high_rewE{c,:,side}];
    kks{7} = [conds.low_rewE{c,:,side}];
    
    indata{1} = out.default{c,side};
    indata{2} = out.fast_e{c,side};
    indata{3} = out.slow_e{c,side};
    indata{4} = out.correct{c,side};
    indata{5} = out.incorrect{c,side};
    indata{6} = out.high_rew{c,side};
    indata{7} = out.low_rew{c,side};
    
    for iii = 1:7
        if ~isempty(kks{iii})
            if size(kks{iii},2)>1
                outdata{iii} = [indata{iii} squeeze(data(chn,:,kks{iii}))];
            else
                outdata{iii} = [indata{iii}];
            end
        else
            outdata{iii} = [indata{iii}];
        end
    end
    
    out.default{c,side} = outdata{1};
    out.fast_e{c,side} = outdata{2};
    out.slow_e{c,side}= outdata{3};
    out.correct{c,side}= outdata{4};
    out.incorrect{c,side}= outdata{5};
    out.high_rew{c,side} = outdata{6};
    out.low_rew{c,side} = outdata{7};
end
end