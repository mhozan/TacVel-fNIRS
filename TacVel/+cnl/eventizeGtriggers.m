function [ALLSTIM] = eventizeGtriggers(triggerAnalog,varargin)
%EVENTIZEGTRIGGERS Automatically tags the onset of every pulse train in the
%analog Galileo Triggers that are recorded by ADI. It also labels the
%onsets with the velocity of each pulsetrain (based on the time interval 
%from fingertips onward). Finally, it creates a Stimulus file readable by
%the Huppert NIRS toolbox.
%   Example(s):
%   NIRSSTIM = eventizeGtriggers(triggerAnalog);
%
% Written by:
% Mohsen Hozan hozan@huskers.unl.edu
% Communication Neuroscience Laboratories
% Center for Brain, Biology, and Behavior
% University of Nebraska-Lincoln
% November 2020
%
% see also
% nirs, Dictionary, nirs.design.StimulusEvents

diff_trig_analog = [0;diff(triggerAnalog(:))];
min_seconds_between_pulsetrains = 20;
minpeak_height= 0.2;

onset_shift = 0;
if nargin >1
    onset_shift = varargin{1};
end




minimum_known_velocity = 4; % cm/s the beginning of each experiment is tagged with a 4cm/s pulse train, 3 minutes before the actual experiment pulse trains (15,20,25,35,45,60) cm/s
% the data recorded before 20201110 does not have the 4cm/s tag, and has the minimum velocity of 20 or 15 cm/s.
%     minvelocities = [...cm/s
%     20; %20200923A_JG
%     20; %20200925A_RW
%     20; %20201002AMH
%     20; %20201021ADM_5ch
%     15;  %20201109AKP
%     4]; %20201111AKP

if nargin == 3
    minimum_known_velocity = varargin{2};
end


[~,onsets] = findpeaks(diff_trig_analog,1000,'MinPeakHeight',minpeak_height);
onsets = onsets+onset_shift;
diff_onsets = [1.1*min_seconds_between_pulsetrains;diff(onsets)];
rising_edge = find(diff_onsets>min_seconds_between_pulsetrains);

trainonsets = onsets(rising_edge);
N_trains = length(trainonsets);

whichvelocity = (onsets(rising_edge+2)-onsets(rising_edge+1));
% whichvelocity = round(2*whichvelocity,2);
% whichvelocity = round(whichvelocity,5);
% [vel_unique_sorted,~,stimvelocityorder] = unique(whichvelocity,'sorted');
[vel_unique_sorted,~,stimvelocityorder] = uniquetol(2*whichvelocity,.01);

N_distinctvelocities = length(vel_unique_sorted);
vel_unique_sorted = vel_unique_sorted./max(vel_unique_sorted);
proportionate_vel_tags = round(minimum_known_velocity./vel_unique_sorted);
aretagsinvalid = ~ismember(proportionate_vel_tags,[4,10:5:70]);
if any(aretagsinvalid)
    warning('Some of the autamically-calculated velocity tags are not a multiple of 5. replacing them with the nearest multiple of 5.')
    proportionate_vel_tags(aretagsinvalid)=round(2*proportionate_vel_tags(aretagsinvalid),-1)/2;    
end
    
velocitylabels = compose("v%02dcmps",proportionate_vel_tags);
sortedtrainlabels = string(stimvelocityorder);

for nn=1:N_trains
    for velcount = 1:N_distinctvelocities
        switch stimvelocityorder(nn)
            case velcount
                sortedtrainlabels(nn) = velocitylabels(velcount);
                st{velcount}.stimonsets(nn) = trainonsets(nn);
        end
    end
end


% taggedvels.onsets = stimonsets;
% taggedvels.labels = stimvelocityorderlabels;


%creating stimulus object for nirs toolbox
keys = cell(1,N_distinctvelocities);
values = cell(1,N_distinctvelocities);
for  ii =1:N_distinctvelocities
    name = char(velocitylabels(ii));
    keys{ii} = name;
    onset = trainonsets(stimvelocityorder==ii);
    dur = 10*ones(size(onset)); %seconds
    amp = 1*ones(size(onset));
    st = nirs.design.StimulusEvents(name,onset,dur,amp);
    
    values{ii} = st;
end

if length(unique(keys))~= length(keys) 
    %THIS SHOULD NOT HAPPEN, 
    %unless two velocities are very close to each other and when rounded, 
    %produce the same integer. if it does happen, you can combine the 2 or
    %give them different tags inside this 'if' statement.
%     minvelocities = [...cm/s
%     20; %20200923A_JG
%     20; %20200925A_RW
%     20; %20201002AMH
%     20; %20201021ADM_5ch
%     15;  %20201109AKP
%     4]; %20201111AKP
    
    s=dbstack;
    warning([mfilename, ' line ', num2str(s(2).line), ': Two different velocities with the same tag.']) 
end

ALLSTIM = Dictionary(keys,values);
