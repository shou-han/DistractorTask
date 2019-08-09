function  addSubjectPlot(group, participant )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
load(group)
allsubjG = allsubj;
CPP_sideG = CPP_side;
CPPr_sideG = CPPr_side;
ERP_sideG = ERPr_side;
N2c_sideG = N2c_side;
N2i_sideG = N2i_side;
RTAG = RTA;
RTsG = RTs;
tG = t;
trG = tr;
hitG = hit;
missG = miss;

load(participant)
allsubj = [allsubjG allsubj];
CPP_side = [CPP_sideG; CPP_side];
CPPr_side = [CPPr_sideG; CPPr_side];
ERPr_side = [ERP_sideG; ERPr_side];
N2c_side = [N2c_sideG; N2c_side];
N2i_side = [N2i_sideG; N2i_side];
RTA.o= [RTAG.o; RTA.o];
RTA.e= [RTAG.e; RTA.e];
RTs = [RTsG;  RTs];
t = [tG t];
tr = [trG tr];
hit = [hitG hit];
miss = [missG miss];

save(group, 'allsubj','CPP_side','CPPr_side','ERPr_side','N2c_side','N2i_side','RTA','RTs','t','tr','chanlocs','side_tags','hit','miss')
end

