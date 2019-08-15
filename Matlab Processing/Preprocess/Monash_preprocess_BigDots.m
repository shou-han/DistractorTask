%% Time

targcodes = [101:148];
nchan = 64;

fs = 500; % new sample rate

%-0.7 seconds before and 5 seconds after
ts = -0.700*fs:1.800*fs; 
t = ts*1000/fs;

BL_time = [-100 0];   % baseline interval in ms
default_response_time = 1.700;

ERP_samps = length(ts);

%% Filters

LPFcutoff_35Hz=35;       % Low Pass Filter cutoff
LPFcutoff_8Hz=8;       % Low Pass Filter cutoff

HPFcutoff=0.1;       % High Pass Filter cutoff

LPF = 1;    % 1 = low-pass filter the data, 0=don't.
HPF = 0;

%% Artifact rejection

ARchans = [1:64];  % just the ones we care about, and RH opposite LH to be symmetric
artifth = 100;
artifchans=[];  % keep track of channels on which the threshold is exceeded, causing trial rejection

chanlocs = readlocs('actiCAP64_2_ThetaPhi.elp','filetype','besa'); %DN for actiCAP
chanlocs = chanlocs';

%% Initialise
% # % Make sure all initialised
numtr=0;
allRT=[]; allrespLR=[]; allTrig=[]; allblock_count = [];
erp_LPF_8Hz = []; erp_LPF_35Hz = []; erp_LPF_8Hz_CSD = []; erp_LPF_35Hz_CSD = [];
% artifacts to reject
artifchans_pretarg = []; artifchans_BL_resp = []; artifchans_resp = []; artifchans_1000ms = [];
pretarg_artrej = []; BL_resp_artrej = []; resp_artrej = []; t1000ms_artrej = [];
ET_pretarg_artrej = []; ET_BL_resp_artrej = []; ET_resp_artrej = []; ET_t1000ms_artrej = [];
all_speed = [];
%% Begin loop
disp(subject_folder{sT})
for f=1:length(files)
    numtrn = 0;
    disp(f)
    EEG = pop_loadbv(paths{f},files{f}); %loads the EEG
    loadbvSK_DN % this takes 64 channels

    EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
    EEG.data = double(EEG.data);
    while 1
        if EEG.event(1).type==1 | EEG.event(1).type==-88
            EEG.event(1) = [];
        else
            break;
        end
    end
    
    load(matfiles{f},'trialCond','par','RTs','conditiondescrip');
    fid = fopen(ET_files{f});
    %% eye tracking
    ET_text = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','Headerlines',22,'ReturnOnError',0);
    fclose(fid);
    for i = 1:size(ET_text{1,3},1)
        if strcmp('GAZE_COORDS',ET_text{1,3}(i))
            screen_res(1) = str2num(cell2mat(ET_text{1,6}(i)))+1
            screen_res(2) = str2num(cell2mat(ET_text{1,7}(i)))+1
            continue
        end
    end
    if screen_res(1)==1024, ranger = 250; elseif screen_res(1)==1280, ranger = 250; else disp(screen_res), keyboard, end
    
    middle = screen_res(1)/2;
    
    if ~exist(ET_matfiles{f}, 'file') %DN: if ET matfile NOT been saved
        FixEyelinkMessages %then calculate and save it now
    end
%%
    trialCond = trialCond+100;

    numev = length(EEG.event);
    % Fish out the event triggers and times
    clear trigs stimes
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end

    % interpolate bad channels
    if ~isempty(badchans)
        EEG.chanlocs = chanlocs;
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end
    %% synchronise the eye tracking with the EEG
    % ET stuff
    load(ET_matfiles{f}) %DN: load the ET mat file
    
    % fix missing trig issue
    EEG_trigs=[]; ET_trigs=[];
    for i = 1:length(EEG.event), EEG_trigs(i) = EEG.event(i).type; end
    for i = 1:length(event), ET_trigs(i) = event(i,2); end
    ET_trigs = ET_trigs(find(ET_trigs>100)); EEG_trigs = EEG_trigs(find(EEG_trigs>100));
    
    if length(ET_trigs)>length(EEG_trigs), last_event = ET_trigs(end-1); 
        if ET_trigs(end)==ET_trigs(end-1), event = event(1:end-2,:); save(ET_matfiles{f},'event','-append'); end
    end
    
    plotter = 0;

    % combine EEG and eye tracking together into one file
    EEG_temp = pop_importeyetracker(EEG,ET_matfiles{f},[first_event ...
        last_event],[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},0,1,0,plotter);
    
    ET_data = EEG_temp.data([nchan+1:nchan+4],:); %synchronised eye tracking data
    clear EEG_temp

    
    % Get rid of mains frequency
    EEG = pop_eegfiltnew(EEG,49,51,[],1,0,0,0);
    EEG_LPF_8Hz = EEG; EEG_LPF_35Hz = EEG; clear EEG;
    
    % First LP Filter
    if LPF
        EEG_LPF_8Hz = pop_eegfiltnew(EEG_LPF_8Hz,0,LPFcutoff_8Hz,[]);
        EEG_LPF_35Hz = pop_eegfiltnew(EEG_LPF_35Hz,0,LPFcutoff_35Hz,[]);
    end

    % First HP Filter
    if HPF 
        EEG_LPF_8Hz = pop_eegfiltnew(EEG_LPF_8Hz,HPFcutoff,0,[]); % filter to 0.1Hz 
        EEG_LPF_35Hz = pop_eegfiltnew(EEG_LPF_35Hz,HPFcutoff,0,[]); % filter to 0.1Hz 
        disp('HPF finished')
    end
    
    % average-reference the whole continuous data (safe to do this now after interpolation and filtering):
    EEG_LPF_8Hz.data = EEG_LPF_8Hz.data - repmat(mean(EEG_LPF_8Hz.data(:,:),1),[nchan,1]);
    EEG_LPF_35Hz.data = EEG_LPF_35Hz.data - repmat(mean(EEG_LPF_35Hz.data(:,:),1),[nchan,1]);

    % extract the epoch
    % find when a trial has started
    targtrigs = [];
    for n=1:length(trigs)
        if any(targcodes(:)==trigs(n))
            targtrigs = [targtrigs n];
        end
    end
    
    
    alltrigs{f} = trigs;
    allrts{f}=RTs;
    
    if trigs(targtrigs(end))==trigs(targtrigs(1))
        motion_on = targtrigs(1:end-1); % GL: indices of trigs when motion on. get rid of last trig, it was a repeat
    else
        motion_on = targtrigs;
    end
    
    for n=1:length(motion_on)
        clear ep_LPF_8Hz ep_LPF_35Hz ep_LPF_8Hz_CSD ep_LPF_35Hz_CSD
        
        % find the response time from the EEG
        locktime = stimes(motion_on(n));
        if motion_on(n)<length(trigs)
            if trigs(motion_on(n)+1)==12 || trigs(motion_on(n)+1)==13 %when they press the button
                response_time = stimes(motion_on(n)+1)-locktime; % time in samples from beginning of motion to response.
                response_time = floor(response_time);
            % threshold of response time (don't know if needed but will
            % see)
                if response_time>default_response_time*fs
                    response_time = default_response_time*fs;
                end
            else
                response_time = default_response_time*fs;
            end
        else
            response_time = default_response_time*fs;
        end
        
        try
            % see only the time in between (default is -0.7 to 1.8)
            ep_LPF_8Hz = EEG_LPF_8Hz.data(:,locktime+ts);   
            ep_LPF_35Hz = EEG_LPF_35Hz.data(:,locktime+ts);
            ep_ET = ET_data(:,locktime+ts);
        catch
            disp('EEG ended too soon2')
            numtr = numtr+1;
            allTrig(numtr) = 0;
            targMotion(numtr) = 0;
            allblock_count(numtr) = f;
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
            
            erp_LPF_8Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_LPF_35Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_LPF_8Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_LPF_35Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            ET_trials(:,:,numtr) = zeros(4,ERP_samps);
            
            continue;
        end

        % new trial
        numtr = numtr+1;
        numtrn = numtrn+1;
        % record baseline amplitude for each channel at different frequencies, ONLY FOR ART REJECT
        BLamp = mean(ep_LPF_8Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); 
        ep_LPF_8Hz_BL = ep_LPF_8Hz - repmat(BLamp,[1,length(t)]);
        
        BLamp = mean(ep_LPF_35Hz(:,find(t>=BL_time(1) & t<BL_time(2))),2); 
        ep_LPF_35Hz_BL  = ep_LPF_35Hz - repmat(BLamp,[1,length(t)]);
        
        % check if the person clicked the trial at all
        ep_test = [find(ts<=(response_time+0.1*fs))];
        if isempty(ep_test)
            disp('Empty epoch for art rejection')
            keyboard
        end
        
        % find which CHANNEL has magnitude over 100 (which is why it uses a
        % 2 for the maximum (along the rows)
        artifchans_thistrial_pretarg = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts<=0))),[],2)>artifth));
        artifchans_thistrial_BL_resp = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts>=-0.1*fs & ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial_resp = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial_1000ms = ARchans(find(max(abs(ep_LPF_35Hz_BL(ARchans,find(ts>=-1*fs & ts<=1*fs))),[],2)>artifth));
        
        artifchans_pretarg = [artifchans_pretarg artifchans_thistrial_pretarg];
        artifchans_BL_resp = [artifchans_BL_resp artifchans_thistrial_BL_resp];
        artifchans_resp = [artifchans_resp artifchans_thistrial_resp];
        artifchans_1000ms = [artifchans_1000ms artifchans_thistrial_1000ms];
        
        % 0 = reject, 1 = keep
        if length(artifchans_thistrial_pretarg) > 0, pretarg_artrej(numtr) = 0; else pretarg_artrej(numtr) = 1; end
        if length(artifchans_thistrial_BL_resp) > 0, BL_resp_artrej(numtr) = 0; else BL_resp_artrej(numtr) = 1; end
        if length(artifchans_thistrial_resp) > 0, resp_artrej(numtr) = 0; else resp_artrej(numtr) = 1; end
        if length(artifchans_thistrial_1000ms) > 0, t1000ms_artrej(numtr) = 0; else t1000ms_artrej(numtr) = 1; end
        
        % scres = 1024 x 768: 512, 384 is middle. 3 deg is 76 pixels. Nope!
        artif_ET_pretarg = find(ep_ET(2,find(ts<=0))<middle-ranger | ep_ET(2,find(ts<=0))>middle+ranger);
        artif_ET_BL_resp = find(ep_ET(2,find(ts>=-0.1*fs & ts<=response_time+0.1*fs))<middle-ranger | ep_ET(2,find(ts>=-0.1*fs & ts<=response_time+0.1*fs))>middle+ranger);
        artif_ET_resp = find(ep_ET(2,find(ts<=(response_time+0.1*fs)))<middle-ranger | ep_ET(2,find(ts<=(response_time+0.1*fs)))>middle+ranger);
        artif_ET_1000ms = find(ep_ET(2,find(ts>=-1*fs & ts<=1*fs))<middle-ranger | ep_ET(2,find(ts>=-1*fs & ts<=1*fs))>middle+ranger);
        
        % 0 = reject, 1 = keep
        if length(artif_ET_pretarg) > 0, ET_pretarg_artrej(numtr) = 0; else ET_pretarg_artrej(numtr) = 1; end
        if length(artif_ET_BL_resp) > 0, ET_BL_resp_artrej(numtr) = 0; else ET_BL_resp_artrej(numtr) = 1; end
        if length(artif_ET_resp) > 0, ET_resp_artrej(numtr) = 0; else ET_resp_artrej(numtr) = 1; end
        if length(artif_ET_1000ms) > 0, ET_t1000ms_artrej(numtr) = 0; else ET_t1000ms_artrej(numtr) = 1; end
            
        
        ep_LPF_8Hz = double(ep_LPF_8Hz);
        ep_LPF_35Hz = double(ep_LPF_35Hz);
        
        ep_LPF_8Hz_CSD = CSD(ep_LPF_8Hz,G_CSD,H_CSD);
        ep_LPF_35Hz_CSD = CSD(ep_LPF_35Hz,G_CSD,H_CSD);
        
        erp_LPF_8Hz(:,:,numtr) = ep_LPF_8Hz;
        erp_LPF_35Hz(:,:,numtr) = ep_LPF_35Hz;
        erp_LPF_8Hz_CSD(:,:,numtr) = ep_LPF_8Hz_CSD;
        erp_LPF_35Hz_CSD(:,:,numtr) = ep_LPF_35Hz_CSD;
        ET_trials(:,:,numtr) = ep_ET(:,:);
        
        allTrig(numtr) = trigs(motion_on(n));
        allblock_count(numtr) = f;


        
        try % change this
            if trigs(motion_on(n)+1)==12 %pressed the left button
                allrespLR(numtr) = 1;
                allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
            elseif trigs(motion_on(n)+1)==13% they pressed the right button
                allrespLR(numtr) = 2;
                allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
            elseif trigs(motion_on(n)+1)==5 % they didn't press the button or too late
                allrespLR(numtr) = 3; % no response, to mark it out from artifact trials.
                allRT(numtr) = nan;
            end
        catch
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
        end

        
    end
%         numtrn
%         size(RTs)
    all_speed = [all_speed RTs*fs];
    destr(f) = numtr;        
end

rejected_trials = length(find(BL_resp_artrej==0));
h1 = figure;
hist(artifchans_BL_resp,[1:nchan]); title([allsubj{sT} ': ' num2str(rejected_trials) ' artifacts = ',num2str(round(100*(rejected_trials/length(allRT)))),'%']) % s from runafew
disp([allsubj{sT},' number of trials: ',num2str(length(find(allRT)))])
saveas(h1, ['Figures/chanrej', num2str(subj) '.png'])

[counts,centers] = hist(artifchans_BL_resp,[1:nchan]);
h1 = figure;
topoplot(counts,chanlocs,'plotchans',[1:nchan],'electrodes','numbers');
title(subject_folder{sT})
pause(1)
saveas(h1, ['Figures/caprej', num2str(subj) '.png'])

erp_LPF_8Hz = single(erp_LPF_8Hz);
erp_LPF_35Hz = single(erp_LPF_35Hz);
erp_LPF_8Hz_CSD = single(erp_LPF_8Hz_CSD);
erp_LPF_35Hz_CSD = single(erp_LPF_35Hz_CSD);

% Baseline erp
baseline_erp = mean(erp_LPF_35Hz(:,find(t>=BL_time(1) & t<=BL_time(2)),:),2);
erp_temp = erp_LPF_35Hz-repmat(baseline_erp,[1,size(erp_LPF_35Hz,2),1]); % baseline full erp
    
erp_temp = squeeze(mean(erp_temp(:,:,find(BL_resp_artrej==1)),3));
h1 = figure
plottopo(erp_temp(:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
    min(min(erp_temp(:,:)))  max(max(erp_temp(:,:)))], ...
    'title',['ERP'],'ydir',1)
saveas(h1, ['Figures/capplot_', num2str(subj) '.png'])


save(['../' path_temp subject_folder{sT} '/' allsubj{sT} '_big_dots_erp_new'],'erp_LPF_8Hz','erp_LPF_35Hz','erp_LPF_8Hz_CSD','erp_LPF_35Hz_CSD', ...
    'allRT','allrespLR','allTrig','allblock_count','t','ET_trials', ...
    'artifchans_pretarg','artifchans_BL_resp','artifchans_resp','artifchans_1000ms', ...
    'pretarg_artrej','BL_resp_artrej','resp_artrej','t1000ms_artrej', ...
    'ET_pretarg_artrej','ET_BL_resp_artrej','ET_resp_artrej','ET_t1000ms_artrej','all_speed')% % % % % %% Longer epochs for alpha

return;