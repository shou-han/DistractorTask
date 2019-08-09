clear all
close all
clc
addpath('/Users/Weibinh/Desktop/Research Project/Data')
load('36_1.mat')
trialNoSqueeze=[];

for i=1:2
    for j=1:60
        x=reinf(i).data(j).main_z(:,1);
        %plot(x)
        s=x-x(1); % zeroing the data 
        %storing
        newD{i,j}=s(~isnan(s));
        %figure
        %plot(newD{i,j})
        %averageForce
        percentageMVC(i,j)=(sum(newD{i,j})/length(newD{i,j}));
        averageForce(i,j)=percentageMVC(i,j)*MVC(:,1);
        %areapertime(i,j)=areaD(i,j)/length(newD{i,j});
        mvc1=tr.reinf;
        sfilt = movmean(newD{i,j},10);
        if max(sfilt)>0.01
            initiateL(i,j) = find(sfilt>(0.01) & sfilt<0.015,1,'last');
            initiateF(i,j) = find(sfilt>(0.01) & sfilt<0.015,1,'first');            
        else
            initiateL(i,j)=1;
            initiateF(i,j)=1;
            trialNoSqueeze=[trialNoSqueeze; i,j];
        end
        slopeS=movingslope(sfilt,200,2); %%%variable to change for smoothing of slope
            initiate = initiateF(i,j);
            searchwindow=400;
        
        [peakGradientNew(i,j) peakTimeNew(i,j)]=max(slopeS(initiate:initiate+searchwindow));
        [peakGradientOld(i,j) peakTimeOld(i,j)]=max(slopeS);
        
        peakTimeNew(i,j) = (peakTimeNew(i,j)-1)+initiate;
        peakTimeOld(i,j) = peakTimeOld(i,j);
    end
end


h=figure
plot(peakTimeOld')

%%


keepfigures(h);

for i=2
    for j=60
        x=reinf(i).data(j).main_z(:,1);
        %plot(x)
        s=x-x(1); % zeroing the data 
        %storing
        newD{i,j}=s(~isnan(s));
        %figure
        %plot(newD{i,j})
        %averageForce
        percentageMVC(i,j)=(sum(newD{i,j})/length(newD{i,j}));
        averageForce(i,j)=percentageMVC(i,j)*MVC(:,1);
        %areapertime(i,j)=areaD(i,j)/length(newD{i,j});
        mvc1=tr.reinf;
        sfilt = movmean(newD{i,j},10);
        if max(sfilt)>0.01
            initiateL(i,j) = find(sfilt>(0.01) & sfilt<0.015,1,'last');
            initiateF(i,j) = find(sfilt>(0.01) & sfilt<0.015,1,'first');            
        else
            initiateL(i,j)=1;
            initiateF(i,j)=1;
            trialNoSqueeze=[trialNoSqueeze; i,j];
        end
        slopeS=movingslope(sfilt,200,2); %%%variable to change for smoothing of slope
            initiate = initiateF(i,j);
            searchwindow=400;
        
        [peakGradientNew(i,j) peakTimeNew(i,j)]=max(slopeS(initiate:initiate+searchwindow));
        [peakGradientOld(i,j) peakTimeOld(i,j)]=max(slopeS);
        
        peakTimeNew(i,j) = (peakTimeNew(i,j)-1)+initiate;
        peakTimeOld(i,j) = peakTimeOld(i,j);
    end
end

figure
plot(peakGradientOld','b')
hold on
plot(peakGradientNew','r')


figure
plot(peakTimeOld')
v =(averageForce)';
m =(mvc1)';
p =(percentageMVC)';
maxF = MVC*308.38;

figure
plot(sfilt)
figure
plot(slopeS)
initiate
hold on