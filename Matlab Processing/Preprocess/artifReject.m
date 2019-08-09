function [artifacts ] = artifReject( structIn,beh,trig_times,eye,epochs,info,artifIn)
%ARTIFREJECT Summary of this function goes here
%   Detailed explanation goes here
% variables
BL_time = epochs.BL_time;
ts = epochs.ts;
t = epochs.t;
numtr = info.numtr;
ARchans = info.ARchans;
artifth = info.artifth;
fs = epochs.fs;
% initialise
artifacts = artifIn;
if ~eye
    pretarg_artrej = artifIn.pretarg_artrej;
    BL_resp_artrej = artifIn.BL_resp_artrej;
    resp_artrej = artifIn.resp_artrej;
    t1000ms_artrej = artifIn.t1000ms_artrej;
    artifchans_pretarg = artifIn.artifchans_pretarg;
    artifchans_BL_resp = artifIn.artifchans_BL_resp;
    artifchans_resp = artifIn.artifchans_resp;
    artifchans_1000ms = artifIn.artifchans_1000ms;
else
    ET_pretarg_artrej = artifIn.ET_pretarg_artrej;
    ET_BL_resp_artrej = artifIn.ET_BL_resp_artrej;
    ET_resp_artrej = artifIn.ET_resp_artrej;
    ET_t1000ms_artrej = artifIn.ET_t1000ms_artrej;
end
%

for n=1:length(trig_times.motion_on)
    numtr = numtr+1;
    Datain = squeeze(structIn.data(:,:,numtr));
    response_time =  beh.response_time(numtr);
    if ~eye
        % record baseline amplitude for each channel at different frequencies, ONLY FOR ART REJECT
        BLamp = mean(Datain(:,find(t>=BL_time(1) & t<BL_time(2))),2);
        Datain_BL = Datain - repmat(BLamp,[1,length(t)]);
        
        
        % find which CHANNEL has magnitude over 100 (which is why it uses a
        % 2 for the maximum (along the rows)
        artifchans_thistrial_pretarg = ARchans(find(max(abs(Datain_BL(ARchans,find(ts<=0))),[],2)>artifth));       
        artifchans_thistrial_BL_resp = ARchans(find(max(abs(Datain_BL(ARchans,find(ts>=-0.1*fs & ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial_resp = ARchans(find(max(abs(Datain_BL(ARchans,find(ts<=(response_time+0.1*fs)))),[],2)>artifth));
        artifchans_thistrial_1000ms = ARchans(find(max(abs(Datain_BL(ARchans,find(ts>=-1*fs & ts<=1*fs))),[],2)>artifth));
        
        % store the artifact channels
        artifchans_pretarg{numtr} = [artifchans_thistrial_pretarg];
        artifchans_BL_resp{numtr} = [artifchans_thistrial_BL_resp];
        artifchans_resp{numtr} = [artifchans_thistrial_resp];
        artifchans_1000ms{numtr} = [artifchans_thistrial_1000ms];
        
        % 0 = reject, 1 = keep
        if length(artifchans_thistrial_pretarg) > 0, pretarg_artrej(numtr) = 0; else pretarg_artrej(numtr) = 1; end
        if length(artifchans_thistrial_BL_resp) > 0, BL_resp_artrej(numtr) = 0; else BL_resp_artrej(numtr) = 1; end
        if length(artifchans_thistrial_resp) > 0, resp_artrej(numtr) = 0; else resp_artrej(numtr) = 1; end
        if length(artifchans_thistrial_1000ms) > 0, t1000ms_artrej(numtr) = 0; else t1000ms_artrej(numtr) = 1; end
    else
        % scres = 1024 x 768: 512, 384 is middle. 3 deg is 76 pixels. Nope!
        middle = info.middle;
        ranger = info.ranger;
        
        ep_ET=Datain;

        artif_ET_pretarg = find(ep_ET(2,find(ts<=0))<middle-ranger | ep_ET(2,find(ts<=0))>middle+ranger);
        artif_ET_BL_resp = find(ep_ET(2,find(ts>=-0.1*fs & ts<=response_time+0.1*fs))<middle-ranger | ep_ET(2,find(ts>=-0.1*fs & ts<=response_time+0.1*fs))>middle+ranger);
        artif_ET_resp = find(ep_ET(2,find(ts<=(response_time+0.1*fs)))<middle-ranger | ep_ET(2,find(ts<=(response_time+0.1*fs)))>middle+ranger);
        artif_ET_1000ms = find(ep_ET(2,find(ts>=-1*fs & ts<=1*fs))<middle-ranger | ep_ET(2,find(ts>=-1*fs & ts<=1*fs))>middle+ranger);
        
        % 0 = reject, 1 = keep
        if length(artif_ET_pretarg) > 0, ET_pretarg_artrej(numtr) = 0; else ET_pretarg_artrej(numtr) = 1; end
        if length(artif_ET_BL_resp) > 0, ET_BL_resp_artrej(numtr) = 0; else ET_BL_resp_artrej(numtr) = 1; end
        if length(artif_ET_resp) > 0, ET_resp_artrej(numtr) = 0; else ET_resp_artrej(numtr) = 1; end
        if length(artif_ET_1000ms) > 0, ET_t1000ms_artrej(numtr) = 0; else ET_t1000ms_artrej(numtr) = 1; end
        
    end
end
if ~eye
    artifacts.pretarg_artrej = pretarg_artrej;
    artifacts.BL_resp_artrej = BL_resp_artrej;
    artifacts.resp_artrej = resp_artrej;
    artifacts.t1000ms_artrej = t1000ms_artrej;
    artifacts.artifchans_pretarg = artifchans_pretarg;
    artifacts.artifchans_BL_resp = artifchans_BL_resp;
    artifacts.artifchans_resp = artifchans_resp;
    artifacts.artifchans_1000ms = artifchans_1000ms;
else
    artifacts.ET_pretarg_artrej = ET_pretarg_artrej;
    artifacts.ET_BL_resp_artrej = ET_BL_resp_artrej;
    artifacts.ET_resp_artrej = ET_resp_artrej;
    artifacts.ET_t1000ms_artrej = ET_t1000ms_artrej;
end
end

