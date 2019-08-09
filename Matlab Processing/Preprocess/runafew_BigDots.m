function runafew_BigDots(subj)
%subj =13;
clearvars -except subj block
cd ..
Start
FilesSort1
cd Preprocess
subj
novel = 0;

close all
clc 
warning off
%%runafew_BigDotsrunafew_BigDots
% 120 for block 2, need to change to 60
%% CSD

load ../beg_vals_old

allbadchans
duds = []; % LK_07_04_14 complete4d the wrong paradigm so will delete this participant later
single_participants = [subj];

% process for each subject
% process for each condition
if ~novel
    for s=1:length(single_participants)
        sT = single_participants(s);
        clear paths files matfiles ET_files ET_matfiles; k=0;
        badchans = allbadchans{sT};
        blocks = allblocks{sT};
        for n=1:length(blocks)
            k=k+1;
            if ismember(subject_folder{s},TCD_bigdots)
                files{k} = filesT{sT,k};
                matfiles{k} = matfilesT{sT,k};
                ET_files{k}=ET_filesT{sT,k};
                ET_matfiles{k} = ET_matfilesT{sT,k};
            elseif ismember(subject_folder{s},Monash_bigdots)
                files{k} = filesT{sT,k};
                paths{k} = pathsT{sT,k};
                matfiles{k} = matfilesT{sT,k};
                ET_files{k}=ET_filesT{sT,k};
                ET_matfiles{k} = ET_matfilesT{sT,k};
            end
        end

        if ismember(subject_folder{s},TCD_bigdots)
            G_CSD = G_TCD;
            H_CSD = H_TCD;
            TCD_preprocess_BigDots
        elseif ismember(subject_folder{s},Monash_bigdots)
            G_CSD = G_monash;
            H_CSD = H_monash;
            
            %
            %Monash_preprocess_BigDots2 %20 only
            Monash_preprocess_BigDots %others only
        end
    end
end
% end