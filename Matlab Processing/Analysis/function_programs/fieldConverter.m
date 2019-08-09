function [ out ] = fieldConverter(in, fields, flag)
if flag
    % if from original
    for iii = 1:length(fields)
        out{iii} = eval(['in.',fields{iii}]);
    end
else
    % if to original
    for iii = 1:length(fields)
        out.(fields{iii})=in{iii};
    end    
end
end

