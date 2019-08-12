% single-trial mediation
clear all
close all
clc
addpath(genpath('function_programs'));
addpath(genpath('../MediationToolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/CSDtoolbox/'));
addpath(genpath('/s/cratch/yn70/ShouHan/Distractors/eeglab13_6_5b/'));
addpath(genpath('../BRAVO2-master/'));
%addpath(genpath('../spm5/'));
% define subject
singletrial=0;

if singletrial; msoffset=4; msBoffset=8; else;msoffset=0; msBoffset=0; end; % to omit the ERPs from the single trial stats
chn = 13;
CSD = 1;
CPPslopeDA=[];CPPslopeDP=[];CPPtrajDA=[];CPPtrajDP=[];
CPPrslopeDA=[];CPPrslopeDP=[];CPPrtrajDA=[]; CPPrtrajDP =[];
CPPtraj=[];CPPrtraj=[];
condsAll = [];
mediationsAll=[];mediationsAll1=[];mediationsAll2=[];
alphaPower = [];
oldnew={'','old'};
subjIndx = 0;
subjCount=0;
for onew=2
    if onew==1
        allsubj = 1:20; ch_CPP=13; subjIndxStart=1;%1:20
    else
        allsubj = [1:18 20:21]; ch_CPP = 13; subjIndxStart=1;
    end
end
for s=1:size(allsubj,2)
    subjCount=subjCount+1;
    sT = allsubj(s)
    if onew==1
        subjIndx = sT;
    else
        subjIndx = sT;
    end
    if CSD
        load(['Data/ERPs/group_plots_erp_diff_CSD_1binsStats_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','mediationsA','chanlocs','t','tr')
        load(['Data/ERPs/group_plots_erp_diff_CSD_1binsStats_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediationsSB','mediationsGB','STFT_time','STFT_timer','mediationsSA','mediationBeta','mediationsG')
    else
        load(['Data/ERPs/group_plots_erp_diff_1binsStats_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','mediationsA','chanlocs','t','tr')
        load(['Data/ERPs/group_plots_erp_diff_1binsStats_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediationsSB','mediationsGB','STFT_time','STFT_timer','mediationsSA','mediationBeta','mediationsG')
    end
    clear ms
    %fields = {'ERP','ERPr','ERPdiff','ERPrdiff'};
    fields = {'ERPdiff','ERPrdiff'};
    mediationsG = rmfield(mediationsG,fields);
    mediationsAll = cat(3,mediationsAll,mediationsS.ERP);
    mediations = mediations;
    alphaPower = [alphaPower; mediationsSA.Alphapower'];
    mediations.Alphapower=[mean(mediationsSA.Alphapower(mediationsSA.c==1)), mean(mediationsSA.Alphapower(mediationsSA.c==2))];
    mediationBeta=mediationBeta;
    ms = fieldnames(mediations);
    msB = fieldnames(mediationBeta);
    msA = fieldnames(mediationsA);
    msG = fieldnames(mediationsG);
    msGB = fieldnames(mediationsGB);
    
    if subjIndx==subjIndxStart
        for i=1:length(msA)
            mediationAllA.(msA{i})=[];
            mediationAllA.subjectAll=[];
        end
        for i =1:length(ms)-msoffset
            mediationAllERP.subjectAll = [];
            mediationAll.subjectAll=[];
            mediationAllERP.(ms{i}) =[];
            mediationAll.(ms{i})=[];
        end
        %combine beta and erp
        for i=4:length(msB)-msBoffset
            mediationAll.(msB{i})=[];
        end
        for i = 1:length(msG)
            mediationAllG.(msG{i})=[];
            mediationAllG.subjectAll=[];
        end
        for i = 1:length(msGB)
            mediationAllGB.(msGB{i})=[];
            mediationAllGB.subjectAll=[];
        end
        for i =1:length(msB)
            mediationABeta.(msB{i})=[];
            mediationABeta.subjectAll=[];
        end
    end
    for i=1:length(msA)
        clear temp
        temp = squeeze(mediationsA.(msA{i}));
        subjects=[];
        subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),2)])];
        mediationAllA.(msA{i}) = [mediationAllA.(msA{i}); temp'];
    end
    mediationAllA.subjectAll= [mediationAllA.subjectAll;  subjects'];
    for i = 1:length(ms)-msoffset
        clear temp tempcond tempZ;
        temp  = squeeze(mediations.(ms{i}));
        tempZ=(temp);
        subjects=[];
        subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),2)])];
        mediationAllERP.(ms{i}) = [mediationAllERP.(ms{i}); temp'];
    end
    mediationAllERP.subjectAll = [mediationAllERP.subjectAll;  subjects'];
    for i = 1:length(msG)
        clear temp tempcond tempZ n;
        subjects = [];
        n = length(size(mediationsG.(msG{i})));
        temp = squeeze(mediationsG.(msG{i}));
        if n<3
            temp = (temp');
            subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),1)])];
            n=1;
        else
            temp = (temp);
            subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),1)])];
        end
        mediationAllG.(msG{i})= cat(n, temp, mediationAllG.(msG{i}));
        
    end
    mediationAllG.subjectAll = [mediationAllG.subjectAll;  subjects'];
    for i = 1:length(msGB)
        clear temp tempcond tempZ;
        subjects = [];
        temp = squeeze(mediationsGB.(msGB{i}));
        temp = (temp);
        if length(size(temp))<3
            mediationAllGB.(msGB{i})= [mediationAllGB.(msGB{i});temp'];
            subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),2)])];
        else
            mediationAllGB.(msGB{i})= cat(3,mediationAllGB.(msGB{i}),temp);
            subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),3)])];
        end
    end
    mediationAllGB.subjectAll = [mediationAllGB.subjectAll;  subjects'];
    for i = 1:length(msB)-msBoffset
        clear temp tempcond tempZ;
        temp  = squeeze(mediationBeta.(msB{i}));
        tempZ=(temp);
        subjects=[];
        subjects = [subjects ;subjIndx*ones([1 size(squeeze(temp),2)])];
        mediationABeta.(msB{i}) = [mediationABeta.(msB{i}); temp'];
    end
    mediationABeta.subjectAll = [mediationABeta.subjectAll;  subjects'];
end
disp('loaded')
% combine the betas with the mediations
for i = 4:length(msB)
    
    mediationAll.(msB{i}) = mediationABeta.(msB{i});
end

for i=1:length(ms)-msoffset
    mediationAll.(ms{i}) = mediationAllERP.(ms{i});
    
end
if singletrial
    mediationAll.Alphapower = alphaPower;
end

mediationAll.subjectAll=mediationAllERP.subjectAll;
msA=fieldnames(mediationAll);

%% save the mediationAll to a large matrix
ms= fieldnames(mediationAll);
alldata=[];
kys = 1;
for i=1:length(ms)
    alldata=[alldata mediationAll.(ms{i})];
    for j=1:size(mediationAll.(ms{i}),2)
        if size(mediationAll.(ms{i}),2)>1
            keys{kys-1+j}= [ms{i} '_time_' num2str(j)];
        else
            keys{kys-1+j}= [ms{i}];
        end
    end
    kys = kys+size(mediationAll.(ms{i}),2);
end
%% omit the ERPs from mediationsG
omitfields = {'CPPrslopeTraj','CPPslopeTraj','ERP','ERPr','erpPeakC','erpPeakI','erpLaterpeak','erprSlope'};
mediationGs = rmfield(mediationAllG, omitfields);
omitfields = {'Beta','Betar','BetaI','BetaIr','STFT','STFTr','ERP','ERPr','betarSlope'};
mediationGBs =  rmfield(mediationAllGB, omitfields);
%% save the file to CSV
patht = '../R_stats/Stats/';
if singletrial
    fields = {'CPPslopeTraj','CPPrslopeTraj','Beta','Betar','BetaI','STFT','STFTr','ERP','ERPr','BetaIr'};
    mediationAll = rmfield(mediationAll,fields);
end
%%
%%
if CSD
    struct2csv (mediationAll,[patht 'mediateCSD' oldnew{onew} '.csv']);
    struct2csv (mediationAllA,[patht 'mediateCSDAcc' oldnew{onew} '.csv']);
    struct2csv (mediationGs,[patht 'mediateGs' oldnew{onew} '.csv']);
    struct2csv (mediationGBs,[patht 'mediateGBs' oldnew{onew} '.csv']);
    save(['slopeDataIndividCSD' oldnew{onew} '.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
else
    struct2csv (mediationAll,[patht 'mediatenoCSD' oldnew{onew} '.csv']);
    struct2csv (mediationAllA,[patht 'mediatenoCSDAcc' oldnew{onew} '.csv']);
    struct2csv (mediationGs,[patht 'mediatenoGs' oldnew{onew} '.csv']);
    save(['slopeDataIndividnoCSD' oldnew{onew} '.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
end



return