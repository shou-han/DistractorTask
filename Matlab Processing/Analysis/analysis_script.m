% save all the plots in a group and the plot can choose whichever one it wants
% participant 9, choose a different set of electrodes
function analysis_script(single_participants)
%single_participants=3;
disp('start')
oldnew = {'', 'old'};
old=1 ;
chP=13;
ch_CPP=13;
ch_lr = 45; % left is all contra using erpN
ch_rl = 51; %right is all ipsi using erpN
% these two (new ones) didn't have the cap aligned so the CPP was a bit too
% high
if ismember(single_participants,[5 14]) && ~old
    ch_CPP=[48 13];
end
consWin = 15;
%% Use Current Source Density transformed erp? 1=yes, 0=no
CSD=1;

addpath(genpath('/scratch/yn70/ShouHan/Distractors/CSDtoolbox/'));
addpath(genpath('/scratch/yn70/ShouHan/Distractors/eeglab13_6_5b/'));
addpath('function_programs/');
eeglab

if old
    load ../beg_vals_old % loads the beginning values
else
    load ../beg_vals
end
chanlocs = readlocs('actiCAP64_2_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
TCD_bigdots = {};
path_temp = '../Data/';
skip_step = 10;
%%
side_tags = {'Left','Right'};

for s2 = 1:length(subject_folder)
    for s = 1:length(TCD_bigdots)
        if strcmp(TCD_bigdots{s},subject_folder{s2})
            TCD_index(s) = s2;
        end
    end
    for s = 1:length(Monash_bigdots)
        if strcmp(Monash_bigdots{s},subject_folder{s2})
            Monash_index(s) = s2;
        end
    end
end
%%
duds = [];
%%
if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
end

%% Define channels, having combined Brain Products and Biosemi data
plot_chans = 1:64;
left_hemi = [1 33 34 4 3 37 36 5 38 6 39 7 9 41 8 40 10 42 11 43 12 15 45 14 44 46 47 16];
right_hemi = [32 63 62 31 30 60 61 27 59 28 58 29 26 56 25 57 21 55 22 54 23 20 51 19 52 50 49 18];
centre_chans = [35 48 2 13 17 64 24 53];
elec_pairs = [1,32;33,63;34,62;4,31;3,30;37,60;36,61;5,...
    27;38,59;6,28;39,58;7,29;9,26;41,56;8,25;40,57;10,21;...
    42,55;11,22;43,54;12,23;15,20;45,51;14,19;44,52;46,50;...
    47,49;16,18];

tester = zeros(64,1);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','labels','plotchans',1:64);

tester = zeros(64,1);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','numbers','plotchans',1:64);

%% Triggers
% ITI,left/right
numch=64;
rtlim=[0.300 1.500];

%%%%%%%%%%%check the channels%%%%%%%%%%%%%%
ch_beta_ci = [8; 25];
ch_contra = [right_hemi;left_hemi];
ch_ipsi = [left_hemi;right_hemi];
BL_erp = [-100,0];
BL_beta = [-100 0];
BL_alpha = [-100 0];
ch_alpha =  [46 47 48 49 50];
% zscore threshold
z_thresh = 3;

%% patches compatible and incompatible with hand
targcodes(:,:,1) =   [101   103   117   119   133   135;...
    102   104   118   120   134   136;...
    105   107   121   123   137   139;...
    106   108   122   124   138   140];
targcodes(:,:,2) =   [ 110   112   126   128   142   144;...
    109   111   125   127   141   143;...
    113   115   129   131   145   147;...
    114   116   130   132   146   148];
%% response time bins parameter and search parameters
no_of_bins = 1;
cs=[1 2; 3 4];
dPA = [1 ;2];
bin_counter = fliplr(1:no_of_bins);
% plotting parameters
for bin = 1:no_of_bins
    rt_bins_tags{bin} = num2str(bin);
end

%% Start loop
for s=1:length(single_participants)
    sT = allsubj{s};
    %% define the time epochs for beta alpha and ERPs
    fs = 500; %resampled to 500 in the runafew
    
    ts = -0.700*fs:1.800*fs;
    t = ts*1000/fs;
    
    % resp-locked erps
    trs = [-.600*fs:fs*.100];
    tr = trs*1000/fs;
    
    % beta times
    STFT_time=[];
    no_of_cycles = 8;
    % get a broad range, e.g. beta, 20 to 35Hz
    %fs_STFT = [25]; % or if you want a particular SSVEP frequency
    %stftlen_STFT = round((1000/fs_STFT*no_of_cycles)/2);
    % for SSVEP frequency make sure it's EXACTLY a particular number of cycles of the frequency.
    % check freq_temp_STFT to make sure SSVEP frequency falls on the range
    fs_STFT = [11:21];%[20:26];
    stftlen_STFT = round((1000/round(median(fs_STFT))*no_of_cycles)/2);
    skip_step = 10;
    cc=1;
    for tt = 1:skip_step:length(ts)-(stftlen_STFT)
        tf = tt:tt+stftlen_STFT-1;
        nfft = length(tf);
        freq_temp_STFT = (0:ceil((nfft+1)/2)-1)*fs/nfft;
        STFT_time(cc) = mean(t(tf));
        cc=cc+1;
    end
    
    % alpha times
    STFT_time_alpha=[];
    % get a broad range, e.g. beta, 20 to 35Hz
    %fs_STFT = [25]; % or if you want a particular SSVEP frequency
    %stftlen_STFT = round((1000/fs_STFT*no_of_cycles)/2);
    % for SSVEP frequency make sure it's EXACTLY a particular number of cycles of the frequency.
    % check freq_temp_STFT to make sure SSVEP frequency falls on the range
    fs_STFT_alpha = [9:11];%[20:26];
    stftlen_STFT_alpha = round((1000/round(median(fs_STFT_alpha))*no_of_cycles)/2);
    skip_step = 10;
    cc=1;
    for tt = 1:skip_step:length(ts)-(stftlen_STFT_alpha)
        tf = tt:tt+stftlen_STFT_alpha-1;
        nfft = length(tf);
        freq_temp_STFT = (0:ceil((nfft+1)/2)-1)*fs/nfft;
        STFT_time_alpha(cc) = mean(t(tf));
        cc=cc+1;
    end
    %% Load the trials
    clear indx1 indx2
    
    load(['../' path_temp subject_folder{s} '/' allsubj{s} '_big_dots_erp_new'],'erp_LPF_8Hz','erp_LPF_8Hz_CSD','erp_LPF_35Hz','erp_LPF_35Hz_CSD','allRT','allrespLR','allTrig','allblock_count',...
        'BL_resp_artrej','ET_BL_resp_artrej','all_speed');
    
    
    if CSD
        erp_all=double(erp_LPF_35Hz_CSD); %20Hz
        erp_beta=double(erp_LPF_35Hz_CSD); %35Hz since looking at beta
    else
        erp_all=double(erp_LPF_35Hz);
        erp_beta=double(erp_LPF_35Hz); %35Hz since looking at beta
    end
    trig_all = allTrig;
    
    ET_BL_resp_artrej_all = ET_BL_resp_artrej;
    BL_resp_artrej_all = BL_resp_artrej;
    
    allRT_all=allRT;
    % Baseline erp
    
    % find accuracy all the trials which are correct and incorrect
    
    indx1(1,:) = find(ismember(allTrig,100+(1:2:48)));
    indx2(1,:) = find(ismember(allTrig,100+(2:2:48)));
    
    allacc_all = allrespLR;
    allacc_all(indx1(allrespLR(indx1)==1))=1;
    allacc_all(indx1(allrespLR(indx1)==2))=0;
    allacc_all(indx1(allrespLR(indx1)==3))=nan;
    
    allacc_all(indx2(allrespLR(indx2)==2))=1;
    allacc_all(indx2(allrespLR(indx2)==1))=0;
    allacc_all(indx2(allrespLR(indx2)==3))=nan;
    
    allacc = allacc_all;
    
    
    %% calculate erps and betas
    % Baseline erp
    baseline_erp = mean(erp_beta(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
    erp_beta = erp_beta-repmat(baseline_erp,[1,size(erp_beta,2),1]); % baseline full erp
    disp('Calculating STFT...')
    erp = double(erp_beta); % chan x time x trial
    STFT = [];
    for trial = 1:size(erp,3)
        cc=1;
        for tt = 1:skip_step:size(erp,2)-(stftlen_STFT)
            tf = tt:tt+stftlen_STFT-1; % define time window
            ep = squeeze(erp(:,tf,trial)); % chop out chan x time window
            nfft = size(ep,2);
            ep = detrend(ep')'; % detrend
            fftx = abs(fft(ep,[],2))./(stftlen_STFT/2);
            fftx = fftx(:,1:ceil((nfft+1)/2));
            ind = find(freq_temp_STFT>fs_STFT(1) & freq_temp_STFT<fs_STFT(end) & freq_temp_STFT~=25);%exclude 25Hz since SSVEP
            %         [~,ind] = min(abs(freq_temp_STFT-fs_STFT)); % if you want SSVEP at one particular frequency
            STFT(:,cc,trial) = mean(fftx(:,ind),2);
            cc=cc+1;
        end
        % asymmetry: right minus left. more positive = more right hemi alph
        beta_asym(right_hemi,:,trial) = (STFT(right_hemi,:,trial)-STFT(left_hemi,:,trial))./...
            ((STFT(right_hemi,:,trial)+STFT(left_hemi,:,trial))/2);
    end
    % Baseline beta
    baseline_beta = mean(STFT(:,find(STFT_time>BL_beta(1) & STFT_time<=BL_beta(2)),:),2);
    beta_TSE_base = STFT-repmat(baseline_beta,[1,size(STFT,2),1]); % baseline full erp
    %make new trials
    for trial = 1:size(erp,3)
        % check if left or right hand, if it's odd, then it's right, else,
        % it is left
        if rem(allTrig(trial),2)
            c=1; % left
        else
            c=2; % right
        end
        beta_CI(left_hemi,:,trial) = beta_TSE_base(ch_contra(c,:),:,trial);
        beta_CI(right_hemi,:,trial) = beta_TSE_base(ch_ipsi(c,:),:,trial);
        beta_CI(centre_chans,:,trial) = beta_TSE_base(centre_chans,:,trial);
        
        beta_hemi_diff_temp(left_hemi,:,trial) = (beta_TSE_base(ch_contra(c,:),:,trial) -  beta_TSE_base(ch_ipsi(c,:),:,trial))./...
            ((beta_TSE_base(ch_contra(c,:),:,trial) + beta_TSE_base(ch_ipsi(c,:),:,trial))/2);
        beta_hemi_diff_temp(right_hemi,:,trial) = (beta_TSE_base(ch_ipsi(c,:),:,trial) -  beta_TSE_base(ch_contra(c,:),:,trial))./...
            ((beta_TSE_base(ch_contra(c,:),:,trial) + beta_TSE_base(ch_ipsi(c,:),:,trial))/2);
        beta_hemi_diff_temp(centre_chans,:,trial) = (beta_TSE_base(centre_chans,:,trial) -  beta_TSE_base(centre_chans,:,trial))./...
            ((beta_TSE_base(centre_chans,:,trial) + beta_TSE_base(centre_chans,:,trial))/2);
    end
    beta_CI_base = beta_CI;
    beta_hemi_diff=beta_hemi_diff_temp;
    
    % alphas
    disp('Calculating pre-target Alpha...')
    alpha_bandlimits = [6 11]; % defining the filter for alpha bandpass.
    [H,G]=butter(4,[2*(alpha_bandlimits(1)/fs) 2*(alpha_bandlimits(2)/fs)]); % alpha bandpass for 500Hz
    
    window = 20; % in samples. Time is double this.
    skip_step_alpha = window/2;
    
    % Alpha time
    alpha_t=[]; cca=1;
    for tt = 1:skip_step_alpha:length(t)-window
        alpha_t(:,cca) = mean(t(tt:tt+window-1));
        cca=cca+1;
    end
    % Alpha Spectrotemporal Evolution a la Thut
    alpha_TSE = []; alpha_asym = [];
    % it is left
    if rem(allTrig(trial),2)
        c=1; % left
    else
        c=2; % right
    end
    for trial = 1:size(erp,3)
        % filtering to alpha
        ep_filt = filtfilt(H,G,squeeze(erp(:,:,trial))')';
        % chop off ends and rectify
        ep_filt = abs(ep_filt(:,find(t>=t(1) & t<=t(end))));
        % smooth
        cca=1;
        for tt = 1:skip_step_alpha:size(ep_filt,2)-window
            alpha_TSE(:,cca,trial) = mean(ep_filt(:,tt:tt+window-1),2);
            cca=cca+1;
        end
        alpha_asym(:,:,trial)=0*alpha_TSE(:,:,trial);
        % asymmetry: right minus left. more positive = more right hemi alph
        alpha_asym(right_hemi,:,trial) = (alpha_TSE(ch_contra(c,:),:,trial)-alpha_TSE(ch_ipsi(c,:),:,trial))./...
            ((alpha_TSE(ch_contra(c,:),:,trial)+alpha_TSE(ch_ipsi(c,:),:,trial))/2);
    end
    
    %make erps contra ipsi
    
    for trial=1:size(erp,3)
        if ismember(allTrig(trial),targcodes(:,:,2))
            side=2;
        else
            side=1;
        end
        erpN(left_hemi,:,trial) = erp(ch_contra(side,:),:,trial); %save as new variable since transformation problems
        erpN(right_hemi,:,trial) = erp(ch_ipsi(side,:),:,trial);
        erpN(centre_chans,:,trial) = erp(centre_chans,:,trial);
    end
    clear erp baseline_erp
    erp = erpN;
    disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(find(allTrig)))])
    
    
    % Baseline alpha
    baseline_alpha = mean(alpha_TSE(:,find(alpha_t>BL_alpha(1) & alpha_t<=BL_alpha(2)),:),2);
    alpha_TSE_base = alpha_TSE-repmat(baseline_alpha,[1,size(alpha_TSE,2),1]); % baseline full erp
    
    % Baseline beta
    baseline_beta = mean(STFT(:,STFT_time>BL_beta(1) & STFT_time<=BL_beta(2),:),2);
    beta_TSE_base = STFT-repmat(baseline_beta,[1,size(STFT,2),1]); % baseline full erp
    
    % Baseline erp
    baseline_erp = mean(erp(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
    erp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp
    
    disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(find(allTrig)))])
    
    
    
    %% make response locked
    % for each trial, shift the erp signature to make RT = 0;
    alpha_tr= -600:skip_step_alpha*2:100;
    %Response locked STFT time in samples
    alpha_trs = -0.6/(skip_step_alpha/fs):.100/(skip_step_alpha/fs);
    STFT_timer= -600:skip_step*2:100;
    %Response locked STFT time in samples
    STFT_timers = -0.6/(skip_step/fs):.100/(skip_step/fs);
    
    erpr = zeros(size(erp,1),length(tr),size(erp,3));
    STFTr = zeros(size(erp,1),length(STFT_timer),size(erp,3));
    alphar = zeros(size(erp,1),length(alpha_tr),size(erp,3));
    alphar_asym= zeros(size(erp,1),length(alpha_tr),size(erp,3));
    validrlock = zeros(1,length(allRT)); % length of RTs.
    % for each trial, shift the erp signature to make RT = 0;
    for n=1:length(allRT)
        [blah,RTsamp] = min(abs(t*fs/1000-allRT(n)));
        [blah,RTsampBeta] = min(abs(STFT_time*fs/1000-allRT(n))); % get the sample point of the RT.
        [blah,RTsampAlpha] = min(abs(alpha_t*fs/1000-allRT(n)));
        if      RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t) &...
                RTsampBeta+STFT_timers(1) >0 & RTsampBeta+STFT_timers(end)<=length(STFT_time) & allRT(n)>0 && ...
                RTsampAlpha+alpha_trs(1) >0 & RTsampAlpha+alpha_trs(end)<=length(alpha_t)  % is the RT larger than 1st stim RT point, smaller than last RT point.
            
            erpr(:,:,n) = erp(:,RTsamp+trs,n);
            STFTr(:,:,n) = beta_TSE_base(:,RTsampBeta+STFT_timers,n);
            betar_hemi_diff(:,:,n) = beta_hemi_diff(:,RTsampBeta+STFT_timers,n); %conditions already includes validrlock
            STFTCIr(:,:,n) = beta_CI_base(:,RTsampBeta+STFT_timers,n);
            alphar(:,:,n) = alpha_TSE_base(:,RTsampAlpha+alpha_trs,n);
            alphar_asym(:,:,n) = alpha_asym(:,RTsampAlpha+alpha_trs,n);
            
            validrlock(n)=1;
        end
    end
    %% define search windows
    clear beta_asym_bins RT_bins RT_bins_all
    
    %perform group trials mediation analysis
    j=0;
    if CSD
        window_searchi = [150,450];
        window_searchc = [150 400];
        window_searchpc = [200 350];
        window_searchall = [700 1200];
        window_slope = [-450 -50];
        window_slopeEarly = [-600 -400];
        window_slopeLate = [-250 -50];
        window_slopeStimEarly = [0 250];
        window_slopeStimLate = [600 800];
        window_betaCOnset = [100,300];
        window_betaCPeak = [400 1200];
        window_betaCSlope = [200 700];
        window_betaRSlope = [-500 -100];
        window_diffSlope=[-350 -50];
        window_betaCclick =2; %this is the number of time increments
        wind=[-0.1 0.1];
        windB = [-0.02, 0.02];
    else
        window_searchi = [150,450];
        window_searchc = [150,400];
        window_searchpc = [150 400];
        window_searchall =  [700 1200];
        window_slope = [-450 -50];
        window_slopeEarly = [-550 -300];
        window_slopeLate = [-300 -50];
        window_betaCOnset = [50,800];
        window_betaCPeak = [400 1200];
        window_betaCSlope = [200 700];
        window_betaRSlope = [-600 -175];
        window_diffSlope=[-400 -100];
        window_betaCclick =2; %this is the number of time increments
        wind=[-0.1,0.1];
        windB = [-0.01, 0.01];
    end
    %% calculate the conditions
    bins_count = zeros([1 no_of_bins]);
    
    for d = 1:size(cs,1)
        for side=1:2
            c=cs(d,:);
            ET_BL_resp_artrej = ET_BL_resp_artrej_all;
            BL_resp_artrej = BL_resp_artrej_all;
            allRT = allRT_all;
            allacc = allacc_all;
            
            c_ind =  find(ismember(allTrig,squeeze(targcodes(c,:,side))));
            
            disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(c_ind))])
            
            allstuff.RTs{s,d} = allRT(c_ind)*1000/fs;
            allstuff.acc{s,d} = allacc(c_ind);
            
            %prepare for speed bins
            bins_count = zeros([1 no_of_bins]);
            for side = 1:2
                for iti = 1:size(targcodes,2)
                    % calcs the indices of the triggers for each appropriate trial type.
                    conds{s,d,iti,side} = find(ismember(allTrig,targcodes(c,iti,side)) & allacc==1 &...
                        allRT>rtlim(1)*fs & allRT<rtlim(2)*fs & BL_resp_artrej==1  & validrlock & ET_BL_resp_artrej ); %
                    
                    %condition for accuracy since need to see it's okay first
                    condsACC{s,d,iti,side} = find(ismember(allTrig,targcodes(c,iti,side)) & ET_BL_resp_artrej ); %                    
                    condsT{s,d,iti,side} = find(ismember(allTrig,targcodes(c,iti,side)) );
                    RTs{s,d,iti,side} = allRT([conds{s,d,iti,side}])*1000/fs;
                end
            end
        end
    end
    % bin within the trials
    for d=1:size(cs,1)
        for side=1:2
            for iti=1:size(targcodes,2)
                RT_temp = [RTs{s,d,iti,side}];
                cond_temp=[conds{s,d,iti,side}];
                counters = ones(1,no_of_bins);
                
                %bin it
                [RT_temp,indx] = sort(RT_temp);
                conds_bin_temp  = cond_temp(indx);
                group_size = floor(length(RT_temp)/no_of_bins);
                rem_size = rem(length(RT_temp),no_of_bins);
                
                binch = find(bins_count==min(bins_count),1);
                for bin =  [1:binch-1]
                    conds_all= conds_bin_temp((bin-1)*group_size+1:bin*group_size);
                    conds_bin{s,d,iti,side,bin} = conds_all;
                    bins_count(bin) = bins_count(bin)+group_size;
                end
                conds_all= conds_bin_temp((binch-1)*group_size+1:(binch)*group_size+rem_size);
                conds_bin{s,d,iti,side,binch} = conds_all;
                bins_count(binch) = bins_count(binch)+group_size+rem_size;
                
                for bin = [binch+1:no_of_bins]
                    conds_all= conds_bin_temp((bin-1)*group_size+rem_size+1:bin*group_size+rem_size);
                    conds_bin{s,d,iti,side,bin} = conds_all;
                    bins_count(bin) = bins_count(bin)+group_size;
                end
            end
        end
    end
    disp(['Subject ',allsubj{s},' Condition Valid Trials: ',num2str(length([conds_bin{s,d,:,:,:}])), ...
        ' = ',num2str(round(100*(length([conds_bin{s,d,:,:,:}]))/length([condsT{s,d,:,:}]))),'%'])
    clear beta_asym_bins RT_bins RT_bins_all
    
    %% create the data for data analysis and plotting
    for c=1:size(cs,1)
        for bin=1:no_of_bins
            j=j+1;
            kk=1;
            % ERPS
            clear coef CPP_smooth CPPr_smooth pretemp N2c_smooth N2i_smooth CPPslope_mean_temp
            % determined from all the N2cs and N2is
            pretemp = find(t==0)-1;
            CPP_mean_temp=[];N2c_mean_temp=[];N2i_mean_temp=[];CPPr_mean_temp=[];RT_temp=[];erp_mean_temp=[];erpr_mean_temp=[];
            N2pc_mean_temp=[];CPPslope_mean_temp=[];
            for side=1:2
                CPP_mean_temp = cat(2, CPP_mean_temp, squeeze(mean(erp(ch_CPP,:,[conds_bin{s,c,:,side,bin}]),1)));
                CPPr_mean_temp = cat(2, CPPr_mean_temp, squeeze(mean(erpr(ch_CPP,:,[conds_bin{s,c,:,side,bin}]),1)));
                for k = [conds_bin{s,c,:,side,bin}]
                    % find slope average
                    CPPslope_mean_temp (:,kk) = movingslope(mean(erp(ch_CPP,:,k),1),200,4);
                    CPPrslope_mean_temp (:,kk) = movingslope(mean(erpr(ch_CPP,:,k),1),200,4);
                    kk=kk+1;
                end
                N2c_mean_temp = squeeze(erp(ch_lr,:,[conds_bin{s,c,:,side,bin}]));
                N2i_mean_temp = squeeze(erp(ch_rl,:,[conds_bin{s,c,:,side,bin}]));
                N2pc_mean_temp= squeeze(erp(ch_lr,:,[conds_bin{s,c,:,side,bin}])-erp(ch_rl,:,[conds_bin{s,c,:,side,bin}]));
                RT_temp = allRT([conds_bin{s,c,:,side,bin}]);
                erp_mean_temp = squeeze(erp(:,:,[conds_bin{s,c,:,side,bin}]));
                erpr_mean_temp = squeeze(erpr(:,:,[conds_bin{s,c,:,side,bin}]));
                
                % store the graphs for plotting
                allBins.Rts{s,c,side,bin} = allRT([conds_bin{s,c,:,side,bin}])*1000/fs;
            end
            allBins.ERP_side(s,c,:,:,bin) = squeeze(mean(erp_mean_temp,3));
            allBins.ERPr_side(s,c,:,:,bin) = squeeze(mean(erpr_mean_temp,3));
            
            allBins.CPP_side(s,c,:,bin) = squeeze(mean(CPP_mean_temp,2));
            allBins.CPPr_side(s,c,:,bin) = squeeze(mean(CPPr_mean_temp,2));
            allBins.CPPslopes_side(s,c,:,bin) = squeeze(mean(CPPslope_mean_temp,2));
            allBins.CPPrslopes_side(s,c,:,bin) = squeeze(mean(CPPrslope_mean_temp,2));
            
            allBins.N2c_side(s,c,:,bin) = squeeze(mean(N2c_mean_temp,2));
            allBins.N2i_side(s,c,:,bin) = squeeze(mean(N2i_mean_temp,2));
            allBins.N2pc_side(s,c,:,bin) = squeeze(mean(N2pc_mean_temp,2));
            
            % Betas
            beta_r_side(s,c,:,:,bin) = squeeze(mean(STFTCIr(1:numch,:,[conds_bin{s,c,:,:,bin}]),3)); %(s, chan, STFT_timers, targetside)
            beta_side(s,c,:,:,bin) = squeeze(mean(beta_CI_base(1:numch,:,[conds_bin{s,c,:,:,bin}]),3));
            % asymmetry betas
            beta_asym_side(s,c,:,:,bin) = squeeze(mean(beta_asym(size(beta_asym,1),:,[conds_bin{s,c,:,:,bin}]),3));
            
            beta_contra(s,c,:,bin) = squeeze(mean(mean(beta_CI_base(ch_beta_ci(1,:),:,[conds_bin{s,c,:,:,bin}]),3),1));
            betar_contra(s,c,:,bin) = squeeze(mean(mean(STFTCIr(ch_beta_ci(1,:),:,[conds_bin{s,c,:,:,bin}]),3),1));
            
            beta_ipsi(s,c,:,bin) = squeeze(mean(mean(beta_CI_base(ch_beta_ci(2,:),:,[conds_bin{s,c,:,:,bin}]),3),1));
            betar_ipsi(s,c,:,bin) = squeeze(mean(mean(STFTCIr(ch_beta_ci(2,:),:,[conds_bin{s,c,:,:,bin}]),3),1));
            
            for side=1:2
                % collect all trials
                beta_contra_all{s,c,bin,side} = mean(beta_CI_base(ch_beta_ci(1,:),:,[conds_bin{s,c,:,side,bin}]),1);
                betar_contra_all{s,c,bin,side} = mean(STFTCIr(ch_beta_ci(1,:),:,[conds_bin{s,c,:,side,bin}]),1);
                beta_ipsi_all{s,c,bin,side} = mean(beta_CI_base(ch_beta_ci(2,:),:,[conds_bin{s,c,:,side,bin}]),1);
                betar_ipsi_all{s,c,bin,side} = mean(STFTCIr(ch_beta_ci(2,:),:,[conds_bin{s,c,:,side,bin}]),1);
            end
            kk=1;
            % beta slopes
            for k = [conds_bin{s,c,:,:,bin}]
                % find slope average
                betac_slope_temp (:,kk) = movingslope(mean(beta_CI_base(ch_beta_ci(1,:),:,k),1),10,2); % ten since it's 20ms increments
                betarc_slope_temp (:,kk) = movingslope(mean(STFTCIr(ch_beta_ci(1,:),:,k),1),10,2);
                betai_slope_temp (:,kk) = movingslope(mean(beta_CI_base(ch_beta_ci(2,:),:,k),1),10,2); % ten since it's 20ms increments
                betari_slope_temp (:,kk) = movingslope(mean(STFTCIr(ch_beta_ci(2,:),:,k),1),10,2);
                kk=kk+1;
            end
            
            betaslope.contra(s,c,:,bin) = squeeze(mean(betac_slope_temp,2));
            betaslope.contraR(s,c,:,bin) = squeeze(mean(betarc_slope_temp,2));
            betaslope.ipsi(s,c,:,bin) = squeeze(mean(betai_slope_temp,2));
            betaslope.ipsiR(s,c,:,bin)  = squeeze(mean(betari_slope_temp,2));
            
            %
            RT_side{s,c,bin} = allRT([conds_bin{s,c,:,:,bin}])*1000/fs;
            
            RT_st_data{c,bin} = RT_side{s,c,bin};
            
            
            
            alpha_side(s,c,:,:,bin) = squeeze(mean(alpha_TSE(:,:,[conds_bin{s,c,:,:,bin}]),3));
            alpha_asym_side(s,c,:,:,bin) = squeeze(mean(alpha_asym(:,:,[conds_bin{s,c,:,:,bin}]),3));
            
            % statistics on distractor and bins and hemi level
            mediations.c(j)=c;
            mediations.bins(j)=bin;
            mediations.RT(j) = mean(RT_temp);
            [mediations.N2c_onset(j)] = obtainOnset_final_subj(-N2c_mean_temp,...
                window_searchc,5,t,0); %the onset will have the same fields as ERPs
            [mediations.N2i_onset(j)] = obtainOnset_final_subj(-N2i_mean_temp,...
                window_searchi,5,t,0); %the onset will have the same fields as ERPs
            
            [mediations.N2cpeak(j), mediations.N2cLatency(j)] = obtainPeak(-N2c_mean_temp,...
                window_searchc,5,t,0); %N2c peak but positive
            [mediations.N2ipeak(j), mediations.N2iLatency(j)] = obtainPeak(-N2i_mean_temp,...
                window_searchi,5,t,0); %N2i peak but positive
            [mediations.N2pcpeak(j), mediations.N2pcLatency(j)] = obtainPeak(-N2pc_mean_temp,...
                window_searchpc,5,t,0); %N2c peak but positive
            [mediations.N2c_laterpeak(j), mediations.N2c_laterpt(j)] = obtainPeak(-N2c_mean_temp, window_searchall,5,t,0); %N2c later peak
            [mediations.N2i_laterpeak(j), mediations.N2i_laterpt(j)] = obtainPeak(N2i_mean_temp, window_searchall,5,t,0); %N2i later peak it's positive
            [mediations.CPPonset(j)] = obtainOnset_final_subj(CPP_mean_temp,...
                [0,500],2,t,0); %CPP onset
            [mediations.CPPpeakl(j), mediations.CPPPeakTime(j)] = obtainPeak(squeeze(CPP_mean_temp),....
                [0,1000],5,t,0); %CPP peak
            [mediations.CPPlevel(j), mediations.CPPrTime(j)] = obtainPeak(squeeze(CPPr_mean_temp),...
                [-500,100],5,tr,1); %CPPr peak
            [mediations.CPPrslope(j), mediations.CPPronset(j)] = obtainOnset_slope(squeeze(CPPr_mean_temp), window_slope, ...
                -450, tr);
            mediations.Alphapower(j)= squeeze(mean(mean(alpha_side(s,c,:,find(alpha_t>-500 & alpha_t<-50),bin),3),4));
            
            
            % mediationsBeta
            mediationBeta.bin(j)=bin;
            mediationBeta.c(j)=c;
            % determined from all the N2cs and N2is
            beta_mean_ctemp=[];beta_mean_itemp=[];betar_mean_ctemp=[];betar_mean_itemp=[];RT_temp=[];
            
            beta_mean_ctemp = squeeze(mean(beta_CI_base(ch_beta_ci(1,:),:,[conds_bin{s,c,:,:,bin}]),1));
            betar_mean_ctemp = squeeze(mean(STFTCIr(ch_beta_ci(2,:),:,[conds_bin{s,c,:,:,bin}]),1));
            RT_temp = allRT([conds_bin{s,c,:,:,bin}]);
            
            mediationBeta.RT(j) = mean(RT_temp);
            [mediationBeta.BetaOnset(j)] = obtainOnset_final_subj(-beta_mean_ctemp,...
                window_betaCOnset,5,STFT_time,0); %beta onset
            [mediationBeta.BetaLevel(j), mediationBeta.betaTime(j)] = obtainPeak(squeeze(-betar_mean_ctemp),...
                [-500,100],5,STFT_timer,1); %CPPr peak
            [mediationBeta.BetarSlope(j), mediationBeta.Betaronset(j)] = obtainOnset_beta_slope(-betar_mean_ctemp,...
                window_betaRSlope, -450, STFT_timer);
        end
        
    end
    %% create single trial data
    %mediationsAccuracy
    clear N2cTemp N2iTemp CPPrTemp tTemp tTempi trTemp
    k=0;
    for c=1:size(cs,1)
        for side=1:2
            for iti=1:size(targcodes,2)
                if ~isempty(condsACC{s,c,iti,side})
                    for jj=1:size(condsACC{s,c,iti,side},2)
                        k=k+1;
                        kk = condsACC{s,c,iti,side}(jj);
                        if mod(kk,2) %odd is left hand
                            mediationsA.hand(k)=1;
                        elseif mod(kk,2)==0
                            mediationsA.hand(k)=2; %right hand
                        end
                        if ismember(c,dPA(1,:))
                            mediationsA.c(k)=1;
                            mediationsA.cong(k)=c;
                        else
                            mediationsA.c(k)=2;
                            mediationsA.cong(k)=c-2;
                        end
                        mediationsA.side(k)=side;
                        mediationsA.iti(k)=iti;
                        mediationsA.hit(k)=allacc(kk);
                        mediationsA.RT(k)=allRT(kk);
                    end
                end
            end
        end
    end
    %%
    warning off
    k=0;
    clear N2cTemp N2iTemp CPPrTemp tTemp tTempi trTemp
    for c=1:size(cs,1)
        for bin=1:no_of_bins
            for side=1:2
                for iti=1:size(targcodes,2)
                    if ~isempty(conds_bin{s,c,iti,side,bin})
                        %size(conds_bin{s,c,iti,side,bin},2)
                        for jj=1:size(conds_bin{s,c,iti,side,bin},2)
                            k=k+1;
                            kk = conds_bin{s,c,iti,side,bin}(jj);
                            
                            if ismember(c,dPA(1,:))
                                mediationsS.c(k)=1;
                                mediationsS.cong(k)=c;
                            else
                                mediationsS.c(k)=2;
                                mediationsS.cong(k)=c-2;
                            end
                            mediationsS.bins(k)=bin;
                            mediationsS.side(k)=side;
                            mediationsS.iti(k)=iti;
                            mediationsS.hit(k)=allacc(kk);
                            if mod(kk,2) %odd is right hand
                                mediationsS.hand(k)=2;
                            elseif mod(kk,2)==0
                                mediationsS.hand(k)=1; %left hand
                            end
                            
                            mediationsS.RT(k)=allRT(kk);
                            
                            %N2c N2i CPP
                            N2cTemp = squeeze(movmean(erp(ch_lr,:,kk),150));
                            N2iTemp = squeeze(movmean(erp(ch_rl,:,kk),150));
                            CPPTemp = squeeze(movmean(mean(erp(ch_CPP,:,kk),1),150));
                            CPPTempOnset = squeeze(mean(erp(ch_CPP,:,kk),150));
                            
                            N2pcTemp=squeeze(movmean(erp(ch_lr,:,kk)-erp(ch_rl,:,kk),150));
                            
                            
                            [blah,RTsamp] = min(abs(t*fs/1000-allRT(kk)));
                            CPPrTemp = CPPTemp(RTsamp+trs);
                            CPPrA = squeeze(mean(erp(ch_CPP,RTsamp+trs,kk),1));
                            
                            [mediationsS.N2cpeak(k), indxN2c] = min(N2cTemp(t>220 & t<400));
                            [mediationsS.N2ipeak(k), indxN2i] = min(N2iTemp(t>220 & t<500));%200 to accomodate for N2c
                            [mediationsS.N2cLaterpeak(k), indxN2cL]  = max(N2cTemp(t>400 & t<900));
                            [mediationsS.N2iLaterpeak(k), indxN2iL]  = max(N2iTemp(t>400 & t<900));
                            
                            [mediationsS.N2pcpeak(k), indxN2pc] = min(N2pcTemp(t>200 & t<400));
                            [mediationsS.N2pcLaterpeak(k), indxN2pcL] = max(N2pcTemp(t>400 & t<800));
                            
                            tTemp = t(t>220 & t<400); mediationsS.N2cLatency(k) = tTemp(indxN2c);
                            tTemp = t(t>220 & t<500); mediationsS.N2iLatency(k) = tTemp(indxN2i);
                            tTemp = t(t>400 & t<900); mediationsS.N2cLLatency(k) = tTemp(indxN2cL);
                            tTemp = t(t>400 & t<900); mediationsS.N2iLLatency(k) = tTemp(indxN2iL);
                            
                            tTemp = t(t>200 & t<400); mediationsS.N2pcLatency(k) = tTemp(indxN2pc);
                            tTemp = t(t>400 & t<900); mediationsS.N2pcLLatency(k) = tTemp(indxN2pcL);
                            
                            %slopes
                            CPPrslopeTemp = movingslope(CPPrTemp,50,4);
                            coef = fitlm(tr(tr>window_slope(1) & tr<window_slope(2)),CPPrTemp(tr>window_slope(1) & tr<window_slope(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsS.CPPrslope(k)=coef.Coefficients.Estimate(2);
                            coef = fitlm(tr(tr>window_slopeEarly(1) & tr<window_slopeEarly(2)),CPPrTemp(tr>window_slopeEarly(1) & tr<window_slopeEarly(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsS.CPPrslopeEarly(k)=coef.Coefficients.Estimate(2);
                            coef = fitlm(tr(tr>window_slopeLate(1) & tr<window_slopeLate(2)),CPPrTemp(tr>window_slopeLate(1) & tr<window_slopeLate(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsS.CPPrslopeLate(k)=coef.Coefficients.Estimate(2);
                            mediationsS.RT(k) = allRT(kk)*1000/fs;
                            mediationsS.CPPrslopeTraj(:,k) = CPPrslopeTemp;
                            mediationsS.CPPrslopeDiff(k) = mediationsS.CPPrslopeLate(k)-mediationsS.CPPrslopeEarly(k);
                            % find the CPP at RT
                            mediationsS.CPPlevel(k) = CPPTemp(t==mediationsS.RT(k));
                            [mediationsS.CPPamplitude(k), indxCPPL] = max(CPPrTemp(tr>-100 & tr<100));
                            tTemp = tr(tr>-100 & tr<100); mediationsS.CPPtime(k)= tTemp(indxCPPL);
                            
                            clear tindx onsetIndx
                            CPPonsetTemp = (CPPTempOnset(t>100 & t< mediationsS.RT(k)));
                            onsetIndx = find((diff(sign(CPPonsetTemp)))>0, 1); tindx = (t(t>100 & t<mediationsS.RT(k)));
                            if isempty(onsetIndx); mediationsS.CPPonset(k)=NaN; else %CPP never go below 1
                                mediationsS.CPPonset(k)= tindx(onsetIndx);
                            end
                            
                            % stim locked slopes
                            CPPslopeTemp = movingslope(CPPTemp,50,4);
                            windowS_slopeEarly=nanmean(mediations.CPPonset) + window_slopeStimEarly;
                            windowS_slopeLate= nanmean(mediations.CPPonset) + window_slopeStimLate;
                            
                            coef = fitlm(t(t>window_slope(1) & t<window_slope(2)),CPPTemp(t>window_slope(1) & t<window_slope(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsS.CPPslope(k)=coef.Coefficients.Estimate(2);
                            coef = fitlm(t(t>windowS_slopeEarly(1) & t<windowS_slopeEarly(2)),CPPTemp(t>windowS_slopeEarly(1) & t<windowS_slopeEarly(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsS.CPPslopeEarly(k)=coef.Coefficients.Estimate(2);
                            coef = fitlm(t(t>windowS_slopeLate(1) & t<windowS_slopeLate(2)),CPPTemp(t>windowS_slopeLate(1) & t<windowS_slopeLate(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsS.CPPslopeLate(k)=coef.Coefficients.Estimate(2);
                            mediationsS.CPPslopeTraj(:,k) = CPPslopeTemp;
                            mediationsS.CPPslopeDiff(k) = mediationsS.CPPslopeLate(k)-mediationsS.CPPslopeEarly(k);
                            
                            clear slopeMax indx CPPslopeTemp
                            
                            mediationsS.ERP(:,:,k) = squeeze(erp(:,:,kk));
                            mediationsS.ERPr(:,:,k) = squeeze(erpr(:,:,kk));
                            mediationsS.ERPdiff(:,:,k) = 0*squeeze(erp(:,:,kk));
                            mediationsS.ERPrdiff(:,:,k) = 0*squeeze(erpr(:,:,kk));
                            mediationsS.ERPdiff(left_hemi,:,k) = squeeze(erp(left_hemi,:,kk)-erp(right_hemi,:,kk));
                            mediationsS.ERPrdiff(left_hemi,:,k) = squeeze(erpr(left_hemi,:,kk)-erpr(right_hemi,:,kk));
                            
                            %mediationsBeta
                            if ismember(c,dPA(1,:))
                                mediationsSB.c(k)=1;
                                mediationsSB.cong(k)=c;
                            else
                                mediationsSB.c(k)=2;
                                mediationsSB.cong(k)=c-2;
                            end
                            mediationsSB.bins(k)=bin;
                            mediationsSB.side(k)=side;
                            mediationsSB.iti(k)=iti;
                            
                            if mod(kk,2) %odd is right hand
                                mediationsSB.hand(k)=2;
                            elseif mod(kk,2)==0
                                mediationsSB.hand(k)=1; %left hand
                            end
                            
                            BetaTemp = squeeze(movmean(mean(beta_CI_base(ch_beta_ci(1,:),:,kk),1),10));
                            BetaTemp = BetaTemp - repmat(mean(BetaTemp(STFT_time>-100 & STFT_time<0)),1,length(BetaTemp));
                            [blah,RTsamp] = min(abs(STFT_time*fs/1000-allRT(kk))); % get the sample point of the RT.
                            BetarTemp = BetaTemp(RTsamp+STFT_timers);
                            
                            coef = fitlm(STFT_timer(STFT_timer>-300 & STFT_timer<-150),BetarTemp(STFT_timer>-300 & STFT_timer<-150),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsSB.BetarSlope(k)=coef.Coefficients.Estimate(2);
                            mediationsSB.RT(k) = allRT(kk)*1000/fs;
                            % find the Betalevel at RT, negative since it
                            % is going down
                            mediationsSB.BetaLevel(k) = BetaTemp(RTsamp);
                            clear tindx onsetIndx
                            BetaOnsetTemp = -(BetaTemp(STFT_time>0 & STFT_time< mediationsSB.RT(k)));%make positive
                            onsetIndx = find(diff(sign(BetaOnsetTemp))>0,1); tindx = (STFT_time(STFT_time>0 & STFT_time<mediationsSB.RT(k)));
                            if isempty(onsetIndx); mediationsSB.BetaOnset(k)=NaN; else
                                mediationsSB.BetaOnset(k)= tindx(onsetIndx);
                            end
                            % ipsi betas, negative since it
                            % is going down
                            BetaITemp = squeeze(movmean(mean(beta_CI_base(ch_beta_ci(2,:),:,kk),1),10));
                            BetaITemp = BetaITemp - repmat(mean(BetaITemp(STFT_time>-100 & STFT_time<0)),1,length(BetaITemp));
                            BetarITemp = BetaITemp(RTsamp+STFT_timers);
                            
                            coef = fitlm(STFT_timer(STFT_timer>-300 & STFT_timer<-150),BetarITemp(STFT_timer>-300 & STFT_timer<-150),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                            mediationsSB.BetarISlope(k)=coef.Coefficients.Estimate(2);
                            % find the CPP at RT
                            mediationsSB.BetaILevel(k) = BetaITemp(RTsamp);
                            
                            clear tindx onsetIndx
                            BetaIOnsetTemp = -fliplr(BetaITemp(STFT_time>0 & STFT_time< mediationsSB.RT(k)));
                            onsetIndx = find(BetaIOnsetTemp>windB(1) &BetaIOnsetTemp<windB(2),1); tindx = fliplr(STFT_time(STFT_time>0 & STFT_time<mediationsSB.RT(k)));
                            if isempty(onsetIndx); mediationsSB.BetaIOnset(k)=NaN; else
                                mediationsSB.BetaIOnset(k)= tindx(onsetIndx);
                            end
                            
                            mediationsSB.Beta(:,k)= squeeze(mean(beta_CI_base(ch_beta_ci(1,:),:,kk),1));
                            mediationsSB.Betar(:,k)= squeeze(mean(beta_CI_base(ch_beta_ci(1,:),RTsamp+STFT_timers,kk),1));
                            mediationsSB.BetaI(:,k)= squeeze(mean(beta_CI_base(ch_beta_ci(2,:),:,kk),1));
                            mediationsSB.BetaIr(:,k)= squeeze(mean(beta_CI_base(ch_beta_ci(2,:),RTsamp+STFT_timers,kk),1));
                            
                            mediationsSB.STFT(:,k)=STFT_time;
                            mediationsSB.STFTr(:,k)=STFT_timer;
                            
                            mediationsSB.ERP(:,:,k)= beta_CI_base(:,:,kk);
                            mediationsSB.ERPr(:,:,k)= beta_CI_base(:,RTsamp+STFT_timers,kk);
                            
                            
                            %mediationsAlpha
                            if ismember(c,dPA(1,:))
                                mediationsSA.c(k)=1;
                                mediationsSA.cong(k)=c;
                            else
                                mediationsSA.c(k)=2;
                                mediationsSA.cong(k)=c-2;
                            end
                            mediationsSA.bins(k)=bin;
                            mediationsSA.side(k)=side;
                            mediationsSA.iti(k)=iti;
                            
                            if mod(kk,2) %odd is right hand
                                mediationsSA.hand(k)=2;
                            elseif mod(kk,2)==0
                                mediationsSA.hand(k)=1; %left hand
                            end
                            
                            %AlphaTemp = squeeze(movmean(mean(alpha_TSE(ch_alpha,:,kk),1),10));
                            %AlphaTemp = AlphaTemp - repmat(mean(AlphaTemp(alpha_t>-100 & alpha_t<0)),1,length(AlphaTemp));
                            %[blah,RTsamp] = min(abs(alpha_t*fs/1000-allRT(kk))); % get the sample point of the RT.
                            %AlpharTemp = AlphaTemp(RTsamp+alpha_trs);
                            
                            % Extract Alpha Power (Pre and Post)
                            mediationsSA.Alphapower(k)= squeeze(mean(mean(mean(alpha_TSE(ch_alpha,find(alpha_t>-500 & alpha_t<-50),kk),1),2),3));
                            
                            mediationsSA.t(:,k)=alpha_t;
                            mediationsSA.tr(:,k)=alpha_tr;
                            
                            mediationsSA.ERP(:,:,k)= alpha_TSE_base(:,:,kk);
                            mediationsSA.ERPr(:,:,k)= alpha_TSE_base(:,RTsamp+alpha_trs,kk);
                            
                        end
                    end
                end
            end
        end
    end
    %%
    warning on
    %% mediation ger style for erps
    kL=1:length(mediationsS.RT);
    mms = fieldnames(mediationsS);
    
    clear mediationsM;
    %mediationsM= mediationsS;
    ks = 1;%size(targcodes,2)
    for iti=1:1 %zscore all the trials
        condsM = kL;%(mediationsS.iti==iti)
        mediationsM.c(:,ks:ks+length(condsM)-1)=mediationsS.c(:,condsM);
        mediationsM.cong(:,ks:ks+length(condsM)-1)=mediationsS.cong(:,condsM);
        mediationsM.iti(:,ks:ks+length(condsM)-1)=mediationsS.iti(:,condsM);
        mediationsM.hand(:,ks:ks+length(condsM)-1)=mediationsS.hand(:,condsM);
        mediationsM.side(:,ks:ks+length(condsM)-1)=mediationsS.side(:,condsM);
        mediationsM.hit(:,ks:ks+length(condsM)-1)=mediationsS.hit(:,condsM);
        mediationsM.bins(:,ks:ks+length(condsM)-1)=mediationsS.bins(:,condsM);
        for f=7:length(mms)
            if length(size(mediationsS.(mms{f})))<3
                mediationsM.(mms{f})(:,ks:ks+length(condsM)-1) = nanzscore(mediationsS.(mms{f})(:,condsM),0,2);
            else
                mediationsM.(mms{f})(:,:,ks:ks+length(condsM)-1) = nanzscore(mediationsS.(mms{f})(:,:,condsM),0,3);
            end
        end
        ks = ks+length(condsM);
    end
    % sort out the trials zscored but within condition so that it's same
    % number of trials
    kk=1;binA=0;
    [blah, indxA]= sort(mediationsM.RT);
    for c=1:size(cs,1)
        c
        clear indx
        indx = indxA(mediationsM.c(indxA)==c);
        bin=0;binA=0;
        grp = floor(length(indx)/no_of_bins);
        for i=1:grp:length(indx)-rem(length(indx),grp)
            clear condsbinc
            bin= bin+1;binA=binA+1 %binA is an index for all the bins
            %         for f=length(mms)-1:length(mms)
            %             mediationsT.(mms{f}) = mediationsM.(mms{f})(:,:,indx(c,:));
            %         end
            condsbinc = indx(grp*(bin-1)+1:grp*(bin));
            
            mediationsG.bin(kk)=binA;
            mediationsG.c(kk)=mean(mediationsM.c(:,condsbinc),2);
            
            % do it for all electrodes
            warning off
            for f=7:length(mms)
                if length(size(mediationsS.(mms{f})))<3
                    clear tempG tempGC datatemp
                    mediationsG.(mms{f})(:,kk)=nanmean(mediationsM.(mms{f})(:,condsbinc),2);
                else
                    mediationsG.(mms{f})(:,:,kk) = nanmean(mediationsM.(mms{f})(:,:,condsbinc),3);
                end
            end
            %ERP_all(:,:,kk) = nanmean(mediationsT.ERP(:,:,condsbinc{c,bin}),3);
            for ch=1:64
                clear tempERP tempERPr tempERPdiff
                tempERP = nanmean(mediationsM.ERP(ch,:,condsbinc),3);
                tempERPr = nanmean(mediationsM.ERPr(ch,:,condsbinc),3);
                
                [mediationsG.erpPeakC(ch,kk), indxN2c] = min(tempERP(t>220 & t<400));
                [mediationsG.erpPeakI(ch,kk), indxN2i] = min(tempERP(t>220 & t<450));
                [mediationsG.erpLaterpeak(ch,kk), indxN2cL] = min(tempERP(t>600 & t<800));
                
                coef = fitlm(tr(tr>window_slope(1) & tr<window_slope(2)),tempERPr(tr>window_slope(1) & tr<window_slope(2)),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                mediationsG.erprSlope(ch,kk)=coef.Coefficients.Estimate(2);
            end
            kk=kk+1;
        end
    end
    mediaitonsGERP = mediationsG;
    
    %% do mediations Ger style for beta
    clear mediationsG mediationsM mediationsT ERP_all
    warning on
    kL=1:length(mediationsSB.RT);
    mms = fieldnames(mediationsSB);
    %mediationsM= mediationsS;
    ks = 1;
    for iti=1%:size(targcodes,2)
        condsM = kL;%(mediationsSB.iti==iti);
        mediationsM.c(:,ks:ks+length(condsM)-1)=mediationsSB.c(:,condsM);
        mediationsM.iti(:,ks:ks+length(condsM)-1)=mediationsSB.iti(:,condsM);
        mediationsM.cong(:,ks:ks+length(condsM)-1)=mediationsSB.cong(:,condsM);
        mediationsM.bins(:,ks:ks+length(condsM)-1)=mediationsSB.bins(:,condsM);
        mediationsM.hand(:,ks:ks+length(condsM)-1)=mediationsSB.hand(:,condsM);
        mediationsM.side(:,ks:ks+length(condsM)-1)=mediationsSB.side(:,condsM);
        for f=5:length(mms)
            if length(size(mediationsSB.(mms{f})))<3
                mediationsM.(mms{f})(:,ks:ks+length(condsM)-1) = nanzscore(mediationsSB.(mms{f})(:,condsM),0,2);
            else
                mediationsM.(mms{f})(:,:,ks:ks+length(condsM)-1) = nanzscore(mediationsSB.(mms{f})(:,:,condsM),0,3);
            end
        end
        ks = ks+length(condsM);
    end
    % sort out the trials zscored but within condition so that it's same
    % number of trials
    clear indxA
    kk=1;binA=0;bin=0;
    [blah, indxA]= sort(mediationsM.RT);
    for c=1:size(cs,1)
        c
        clear indx
        indx = indxA(mediationsM.c(indxA)==c); bin=0;binA=0;        grp = floor(length(indx));
        for i=1:grp:length(indx)-rem(length(indx),grp)
            clear condsbinc
            bin= bin+1;binA=binA+1; %binA is an index for all the bins
            %         for f=length(mms)-1:length(mms)
            %             mediationsT.(mms{f}) = mediationsM.(mms{f})(:,:,indx(c,:));
            %         end
            condsbinc = indx(grp*(bin-1)+1:grp*(bin));
            
            mediationsG.bin(kk)=binA;
            mediationsG.c(kk)=mean(mediationsM.c(:,condsbinc),2);
            
            % do it for all electrodes
            warning off
            for f=7:length(mms)
                if length(size(mediationsSB.(mms{f})))<3
                    clear tempG tempGC datatemp
                    mediationsG.(mms{f})(:,kk)=nanmean(mediationsM.(mms{f})(:,condsbinc),2);
                else
                    mediationsG.(mms{f})(:,:,kk) = nanmean(mediationsM.(mms{f})(:,:,condsbinc),3);
                end
            end
            %ERP_all(:,:,kk) = nanmean(mediationsT.ERP(:,:,(st(c,bin):ed(c,bin))),3);
            for ch=1:64
                tempERP = nanmean(mediationsM.ERP(ch,:,condsbinc),3);
                tempERPr = nanmean(mediationsM.ERPr(ch,:,condsbinc),3);
                coef = fitlm(STFT_timer(STFT_timer>-300 & STFT_timer<-150),tempERPr(STFT_timer>-300 & STFT_timer<-150),'RobustOpts','on');% coef gives 2 coefficients fitting r = slope * x + intercept
                mediationsG.betarSlope(ch,kk)=coef.Coefficients.Estimate(2);
            end
            kk=kk+1;
        end
    end
    warning on
    mediationsGB = mediationsG;
    mediationsG = mediaitonsGERP;
end

%% save the data for faster plots
if CSD
    save(['Data/ERPs/group_plots_erp_diff_CSD_1binsStats_' oldnew{old+1} '_' num2str(single_participants) '_' num2str(chP)],'RTs','t','plot_chans','tr','chanlocs',...
        'allBins','allsubj','subject_folder', 'allstuff','no_of_bins','mediations','mediationsS','mediationsA','mediationsG','mediationsSA',...
        'conds','condsACC','conds_bin',...
        'beta_side','beta_r_side','STFT_time','STFT_timer',...
        'beta_contra','betar_contra','beta_ipsi','betar_ipsi','RT_side','beta_hemi_diff','betar_hemi_diff',...
        'mediationBeta','mediationsSB','betaslope','mediationsGB','alpha_side','alpha_asym_side','alpha_t','alpha_tr');
else
    save(['Data/ERPs/group_plots_erp_diff_1binsStats_' oldnew{old+1} '_' num2str(single_participants) '_' num2str(chP)],'RTs','t','plot_chans','tr','chanlocs',...
        'allBins','allsubj','subject_folder', 'allstuff','no_of_bins','mediations','mediationsS','mediationsA','mediationsG','mediationsSA',...
        'conds','condsACC','conds_bin',...
        'beta_side','beta_r_side','STFT_time','STFT_timer',...
        'beta_contra','betar_contra','beta_ipsi','betar_ipsi','RT_side','beta_hemi_diff','betar_hemi_diff',...
        'mediationBeta','mediationsSB','betaslope','mediationsGB','alpha_side','alpha_asym_side','alpha_t','alpha_tr');
end


%Plot_original