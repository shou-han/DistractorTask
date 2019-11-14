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
noBins=1;%%%%%%%%%%%%%%%to change
for s=1:size(allsubj,2)
    subjCount=subjCount+1;
    sT = allsubj(s)
    if onew==1
        subjIndx = sT;
    else
        subjIndx = sT;
    end
    if CSD
        load(['Data/ERPs/group_plots_erp_diff_CSD_5binsStats_8Hz_allTrials_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','mediationsA','chanlocs','t','tr')
        load(['Data/ERPs/group_plots_erp_diff_CSD_5binsStats_35Hz_allTrials_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediationsSB','mediationsGB','STFT_time','STFT_timer','mediationsSA','mediationBeta','mediationsG')
    else
        load(['Data/ERPs/group_plots_erp_diff_5binsStats_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediations','mediationsS','mediationsA','chanlocs','t','tr')
        load(['Data/ERPs/group_plots_erp_diff_5binsStats_' oldnew{onew} '_' num2str(sT) '_' num2str(chn)],'mediationsSB','mediationsGB','STFT_time','STFT_timer','mediationsSA','mediationBeta','mediationsG')
    end
    clear ms
    fields2 = {'ERP','ERPr','ERPdiff','ERPrdiff'};
    fields = {'ERPdiff','ERPrdiff'};
    mediationsG = rmfield(mediationsG,fields);
    mediationsAll = cat(3,mediationsAll,mediationsS.ERP);
    mediationsS.RT = zscore(mediationsS.RT);
    mediations = mediationsS;mediations= rmfield(mediations,fields2);
     %size(mediations.CPPonset)
    alphaPower = [alphaPower; mediationsSA.Alphapower'];
%     for c=1:2
%     for bins=1:noBins
%     mediations.Alphapower(bins+(c-1)*noBins)=[mean(mediationsSA.Alphapower(mediationsSA.c==c & mediationsSA.bins==bins))];
%     end
%     end
    mediations.Alphapower=mediationsSA.Alphapower;
    fieldsB = {'ERP','ERPr','STFT','STFTr','BetarISlope','BetaILevel','BetaIOnset','cong','bins','side','RT','Beta','Betar','BetaI','BetaIr'};
    mediationBeta= rmfield(mediationsSB,fieldsB);
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
mediationAll.CPPonset(mediationAll.CPPonset==5000)=NaN;
%% remove unnecessary fields
patht = '../R_stats/Stats/';
if singletrial
    fields = {'CPPslopeTraj','CPPrslopeTraj','Beta','Betar','BetaI','STFT','STFTr','ERP','ERPr','BetaIr'};
    mediationAll = rmfield(mediationAll,fields);
end

%% order the z score and divide them into DADP and then divide them into eight bins, 
no_of_bins=8;
mms = fieldnames(mediationAll);
clear mediationBins;
[blah, indxA]= sort(mediationAll.RT);kk=1;
for c=1:2
    indx = indxA(mediationAll.c(indxA)==c);
    grp = floor(length(indx)/no_of_bins);
    for bin=1:no_of_bins
        condsbinc = indx(grp*(bin-1)+1:grp*(bin));
        subjInclude = unique(mediationAll.subjectAll(condsbinc,:));
        clear condsbincsubj
        for s = 1:length(subjInclude)
            %mediationBins.c(kk,:)=mean(mediationAll.c(condsbinc,:),1);
            subjs=subjInclude(s);
            condsbincsubj = condsbinc(mediationAll.c(condsbinc)==c & mediationAll.subjectAll(condsbinc)==subjs);
            mediationBins.c(kk,:)=mean(mediationAll.c(condsbincsubj,:),1);
            mediationBins.count(kk,:)=length(condsbincsubj);
            % do it for all electrodes
            warning off
            for f=1:length(mms)
                clear tempG tempGC datatemp
                mediationBins.(mms{f})(kk,:)=nanmean(mediationAll.(mms{f})(condsbincsubj,:),1);
            end
            mediationBins.bins(kk,:)=bin;
            kk=kk+1;
        end
    end
end
omitfields = {'CPPrslopeTraj','CPPslopeTraj'};
mediationBinsA = rmfield(mediationBins, omitfields);
for bin=1:no_of_bins
bar(bin,mean(mediationBinsA.CPPrslope(mediationBinsA.c==1 & mediationBinsA.bins==bin)),0.4,'b');
hold on;
bar(bin+0.4,mean(mediationBinsA.CPPrslope(mediationBinsA.c==2 & mediationBinsA.bins==bin)),0.4,'r');
end
figure
for bin=1:no_of_bins
bar(bin,mean(mediationBinsA.count(mediationBinsA.c==1 & mediationBinsA.bins==bin)),0.4,'b');
hold on;
bar(bin+0.4,mean(mediationBinsA.count(mediationBinsA.c==2 & mediationBinsA.bins==bin)),0.4,'r');
end
%% order the z score and divide them into eight bins, and then divide them into DA DP
no_of_bins=8;
mms = fieldnames(mediationAll);
[blah, indxA]= sort(mediationAll.RT);kk=1;
grp = floor(length(indxA)/no_of_bins);
for bin=1:no_of_bins
    indx = indxA(grp*(bin-1)+1:grp*(bin));
    for c=1:2
        condsbinc = indx(mediationAll.c(indx)==c);
        subjInclude = unique(mediationAll.subjectAll(condsbinc,:));
        clear condsbincsubj
        for s = 1:length(subjInclude)
            subjs=subjInclude(s);
            condsbincsubj = condsbinc(mediationAll.c(condsbinc)==c & mediationAll.subjectAll(condsbinc)==subjs);
            mediationBins.c(kk,:)=mean(mediationAll.c(condsbincsubj,:),1);
            mediationBins.count(kk,:)=length(condsbincsubj);
            % do it for all electrodes
            warning off
            for f=1:length(mms)
                clear tempG tempGC datatemp
                mediationBins.(mms{f})(kk,:)=nanmean(mediationAll.(mms{f})(condsbincsubj,:),1);
            end
            mediationBins.bins(kk,:)=bin;
            kk=kk+1;
        end
    end
end
omitfields = {'CPPrslopeTraj','CPPslopeTraj'};
mediationBins = rmfield(mediationBins, omitfields);
figure
for bin=1:no_of_bins
bar(bin,mean(mediationBins.CPPrslope(mediationBins.c==1 & mediationBins.bins==bin)),0.4,'r');
hold on;
bar(bin+0.4,mean(mediationBins.CPPrslope(mediationBins.c==2 & mediationBins.bins==bin)),0.4,'b');
end
figure
for bin=1:no_of_bins
bar(bin,mean(mediationBins.count(mediationBins.c==1 & mediationBins.bins==bin)),0.4,'b');
hold on;
bar(bin+0.4,mean(mediationBins.count(mediationBins.c==2 & mediationBins.bins==bin)),0.4,'r');
end
%%
%%
if CSD
   struct2csv (mediationBins,[patht 'mediateAllZ5binCSD' oldnew{onew} '.csv']);
   struct2csv (mediationAll,[patht 'mediateAllZ5binCSDAcc' oldnew{onew} '.csv']);
   struct2csv (mediationGs,[patht 'mediateAllZ5binGs' oldnew{onew} '.csv']);
   struct2csv (mediationGBs,[patht 'mediateAllZ5binGBs' oldnew{onew} '.csv']);
   save(['slopeDataIndividCSD' oldnew{onew} '.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
else
    struct2csv (mediationAll,[patht 'mediatenoCSD' oldnew{onew} '.csv']);
    struct2csv (mediationAllA,[patht 'mediatenoCSDAcc' oldnew{onew} '.csv']);
    struct2csv (mediationGs,[patht 'mediatenoGs' oldnew{onew} '.csv']);
    save(['slopeDataIndividnoCSD' oldnew{onew} '.mat'],  'CPPslopeDA', 'CPPslopeDP', 'CPPrslopeDA', 'CPPrslopeDP',  't', 'tr', 'mediationAllG', 'mediationAllGB', 'STFT_time', 'STFT_timer')
end



return