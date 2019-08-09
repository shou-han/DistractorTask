function [ ET_out, infoOut] = etload( ETfiles,ET_matfiles,EEG,nchan,info)
%ETLOAD Summary of this function goes here
%   Detailed explanation goes here
    infoOut = info;
    ETfiles
    fid = fopen(ETfiles);
    ET_text = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','Headerlines',22,'ReturnOnError',0);
    fclose(fid);
    for i = 1:size(ET_text{1,3},1)
        if strcmp('GAZE_COORDS',ET_text{1,3}(i))
            screen_res(1) = str2num(cell2mat(ET_text{1,6}(i)))+1
            screen_res(2) = str2num(cell2mat(ET_text{1,7}(i)))+1
            continue
        end
    end
    screen_res = [1024];
    if screen_res(1)==1024, ranger = 76; elseif screen_res(1)==1280, ranger = 98; else disp(screen_res), keyboard, end
    ranger = 76;
    %
    middle = screen_res(1)/2;
    
    if ~exist(ET_matfiles, 'file') %DN: if ET matfile NOT been saved
        FixEyelinkMessages %then calculate and save it now
    end
    
    load(ET_matfiles) %DN: load the ET mat file
    % synchronize the events
    EEG_trigs=[]; ET_trigs=[];
    for i = 1:length(EEG.event), EEG_trigs(i) = EEG.event(i).type; end
    for i = 1:length(event), ET_trigs(i) = event(i,2); end
    ET_trigs = ET_trigs(find(ET_trigs>100)); EEG_trigs = EEG_trigs(find(EEG_trigs>100));
    
    if length(ET_trigs)>length(EEG_trigs), last_event = ET_trigs(end-1); 
        if ET_trigs(end)==ET_trigs(end-1), event = event(1:end-2,:); save(ET_matfiles{f},'event','-append'); end
    end
    
    plotter = 0;
    %Add an extra 4 rows into the a temporary EEG structure
    % - 'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'. 
    % This will add these as extra channels onto %EEG.data. 
    % So the final channel is the pupil area (i.e. diameter):   
    EEG_temp = pop_importeyetracker(EEG,ET_matfiles,[first_event ...
        last_event],[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},0,1,0,plotter);
    [output_cell,~,~] = command_window_text();
    text = output_cell{length(output_cell)-1};
    numsamp = sscanf(text,'%*s%*s%*s%*s%*s%*f%*s%*s%f');
    if numsamp<30
        beep
        disp([allsubj{s},', block ',num2str(f),': ET sync issue'])
        figure, plot(ET_trigs), hold on, plot(EEG_trigs)
        keyboard
    end
    
    ET_data = EEG_temp.data([nchan+1:nchan+4],:);
%     scres = 1280 x 1024: 640,512 is middle
%     temp = ET_data([2,3],[stimes(1):end]); ET_data_mean(1) = mean(temp(1,find(temp(1,:)>0)),2); ET_data_mean(2) = mean(temp(2,find(temp(2,:)>0)),2);
%     ET_data(2,find(ET_data(2,:)>0)) = ET_data(2,find(ET_data(2,:)>0))-ET_data_mean(1)+640;
%     ET_data(3,find(ET_data(3,:)>0)) = ET_data(3,find(ET_data(3,:)>0))-ET_data_mean(2)+512;
    ET_out.data = ET_data;
    infoOut.middle = middle;
    infoOut.ranger = ranger;
end

