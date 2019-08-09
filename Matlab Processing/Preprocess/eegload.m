function [ EEG ] = eegload( paths, files,nchan,chanlocs,badchans)

    EEG = pop_loadbv(paths,files);
    loadbvSK_DN % this takes 65 channels
    EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
    EEG.data = double(EEG.data);
        
    while 1
        if EEG.event(1).type==1 | EEG.event(1).type==-88
            EEG.event(1) = [];
        else
            break;
        end
    end
    
    numev = length(EEG.event);
    % Fish out the event triggers and times
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end

    % interpolate bad channels
    if ~isempty(badchans)
        EEG.chanlocs = chanlocs;
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end    
end

