function [ EEG_p_data ] = eegprocess( EEG,nchan)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %EEGPROCESS 
    LPFcutoff_35Hz=35;       % Low Pass Filter cutoff
    LPFcutoff_8Hz=8;       % Low Pass Filter cutoff

    HPFcutoff=0.1;       % High Pass Filter cutoff

    LPF = 1;    % 1 = low-pass filter the data, 0=don't.
    HPF = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get rid of mains frequency
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
    
    EEG_p_data.Hz8 = EEG_LPF_8Hz;
    EEG_p_data.Hz35 = EEG_LPF_35Hz;
end

