function [ALLSTIM] = eventizeGtriggers2(triggerAnalog,varargin)
%EVENTIZEGTRIGGERS2 Automatically tags the onset of every pulse train in the
%analog Galileo Triggers that are recorded by ADI. It also labels the
%onsets with the velocity of each pulsetrain (based on the time interval 
%from fingertips onward). Finally, it creates a Stimulus file readable by
%the Huppert NIRS toolbox.
% rewritten for the 2021 experiment protocol 20, 25, 31 and 39 cm/s
% The nonlinear separation of the distance between stimulation velocities 
% allows for an equidistance logarithmic spread of the temporal dynamics 
% of the stimuli, in accordance with the Weber'ss law (Nover et al., 2005) 
% when it comes to the perception of sensory stimuli.
% Each recording starts with an equidistant pulstrain lasting 5 seconds,
% followed by at least 90 seconds of silence (baseline recording), 
% then followed by an equidistant dummy trial (15 s on + 25 s  off), 
% which must be exluded from statistical analysis to account for the 
% startle/surprise effect. 
% 40-s Stimulation Trials (15-on 25-off) immediately follow the dummy trial 
% in a randomized order at 4 different velocities of 20, 25, 31, 39 cm/s. 
% The first session ends after 1 dummy + 36 trials (1480s ~= 25 minutes)
% A second identical session is resumed after 5-15 minutes of the
% participant resting, chatting, stretching, and having water/snack.
% Each session is preceded by ~90 s of pre-baseline 
% and followed by ~60s or more of post-baseline,
% i.e., a prolong period of no stimulation which can be used for cleaning
% the data or as a reference period for analyses.
% the recording is uninterrupted throughout, therefore the intermission
% recordings should not be analyzed, except for investigational
% purposes. 
% For some participants the seconds session has the 5s intro, for some it
% doesn't. Therefire the dummy onset + 40s should be used as the timestamp
% for the beginning of the trials. 
% the experiment concludes after the 2nd session. 
% 1st session: 5s intro           + 90s pre-baseline + 1 dummy trials + 36 trials + 60s post-baseline
% 2nd session: 5s intro(optional) + 90s pre-baseline + 1 dummy trials + 36 trials + 60s post-baseline
% 72 total trials, 18 per velocity.

%   Example(s):
%   NIRSSTIM = eventizeGtriggers2(triggerAnalog);
%
% Written by:
% Mohsen Hozan hozan@huskers.unl.edu
% Communication Neuroscience Laboratories
% Center for Brain, Biology, and Behavior
% University of Nebraska-Lincoln
% May 2023
%
% see also
% nirs, Dictionary, nirs.design.StimulusEvents

diff_trig_analog = [0;diff(triggerAnalog(:))];
min_seconds_between_pulsetrains = 20; % each trial is 40 seconds, first 15 s is stim on followed by 25s stim off.
minpeak_height= 0.2;

onset_shift = 0;
if nargin >1
    onset_shift = varargin{1};
end

minimum_known_velocity = 20; % For the 2021 dataset, 20cmps is the lowest velocity.

if nargin >2 
    minimum_known_velocity = varargin{2};
end

Fs = 1000; %1000Hz assumed frequency of triggeranalog, i.e. the auxiliary ADI recordings
if nargin>3
    Fs=varargin{3};
end

trial_onset_shift = 0;  % %stim_on duration is 15 seconds within each 40s trial; this should be 0 by default, but the onset can be optionally adjusted forward or backward.
if nargin>4
    trial_onset_shift=varargin{4};
end

stim_effect_dur = 25; %stim_off duration is 25 seconds within each 40s trial; if the trial_onset_shift is greater than 0, this should be less than 40, so that the sum is never greater than 40 to avoid overlapping with the subsequent trial.
if nargin>5
    stim_effect_dur=varargin{5};
end

trial_dur = stim_effect_dur+trial_onset_shift;
if trial_dur>40
    fprintf(['Trial Duration exceeds 40 seconds. trial_dur=%2.0f seconds.' ...
        '\n Double check the inputs of nirs.cnl.eventizeGtriggers2.\n'],trial_dur)
end


[~,onsets] = findpeaks(diff_trig_analog,Fs,'MinPeakHeight',minpeak_height);
onsets = onsets+onset_shift;
diff_onsets = [min_seconds_between_pulsetrains;diff(onsets)];
rising_edge = find(diff_onsets>min_seconds_between_pulsetrains);

trainonsets = onsets(rising_edge);
N_trains = length(trainonsets);
session_dur = 70+N_trains*40+50; %baseline pre + 37 *40s trials + baseline post ; should work for most participants, unless something has prolonged or interrupted the session.

if N_trains<37
    warning(['%2.0f stimulation trials only\n Is this a 2021 Galileo protocol consisting of 72 trials ' ...
        '(9 stim (x4 velocity per session) + 1 intro?'],N_trains)
end


whichvelocity = (onsets(rising_edge+2)-onsets(rising_edge+1));
% whichvelocity2 = (onsets(rising_edge+3)-onsets(rising_edge+2)); %this gives a bug of 75 dummycount for 20211026BMH
whichvelocity3 = (onsets(rising_edge+4)-onsets(rising_edge+3));
equidistant = abs(whichvelocity3-whichvelocity)<0.01; %equidistant triggers; dummy or intro trial;

% dummycount = length(find(equidistant));
% if dummycount==2
%     warning('check the intro/dummy stimulation')
% elseif dummycount==3
%     %normal
% elseif dummycount>4
%     keyboard
%     % warning('dummycount')
% else
%     %
% end

whichvelocity(equidistant) = 0.03; %a very small arbitrary value to tag the intro and dummy trials

% whichvelocity = round(2*whichvelocity,2);
% whichvelocity = round(whichvelocity,5);
% [vel_unique_sorted,~,stimvelocityorder] = unique(whichvelocity,'sorted');
[vel_unique_sorted,~,stimvelocityorder] = uniquetol(whichvelocity,.02);


tvec = onset_shift + [0:length(triggerAnalog)-1]./Fs;
figure(round(onsets(1)*1000)), clf, plot(tvec,triggerAnalog), hold on, scatter(trainonsets ,1.0+0.1*stimvelocityorder.*ones(size(trainonsets)),'filled')
ylim([-.1,2.5]),drawnow

n_intro = nnz(stimvelocityorder==1);
if n_intro~=1
    warning('%2.0d into stimuli,check if the intro/dummy stimuli are tagged correctly',n_intro)
end

N_distinctvelocities = length(vel_unique_sorted);
vel_unique_sorted = vel_unique_sorted./max(vel_unique_sorted);
proportionate_vel_tags = round(minimum_known_velocity./vel_unique_sorted);
aretagsinvalid = ~ismember(proportionate_vel_tags,[20, 25, 31, 39]);
if sum(aretagsinvalid)>1 %the intro tag is an arbitrary number greater than 100. there should be no other unknown velocities. 
    warning('Unknown Velocity!!!')
    %20221130A session 1 has extra stimuli at the begininng and a lone
    %trigger at 1686.84s
    disp(proportionate_vel_tags(aretagsinvalid).')
    % keyboard
end
    
velocitylabels = compose("v%02dcmps",proportionate_vel_tags);
velocitylabels(strlength(velocitylabels)>7) = "intro";

sortedtriallabels = string(stimvelocityorder);

for nn=1:N_trains
    for velcount = 1:N_distinctvelocities
        switch stimvelocityorder(nn)
            case velcount
                sortedtriallabels(nn) = velocitylabels(velcount);
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
    % dur = (stim_off_dur+stim_on_dur)*ones(size(onset)); %seconds

    onset = onset + trial_onset_shift; %onset after the end of the stimulation
    dur = (stim_effect_dur)*ones(size(onset)); %seconds


    % dur = stim_off_dur*ones(size(onset)); %seconds
    amp = 1*ones(size(onset));
    st = nirs.design.StimulusEvents(name,onset,dur,amp);
    if strcmpi (name,'intro')
        st.regressor_no_interest = true;

        session_onset = onset(1)-70;
        if session_onset<session_dur
            session_name = 'Session1';
        else
            session_name = 'Session2';
        end
        % session_dur = session_dur;
        session_amp = 0.1;
        session_st = nirs.design.StimulusEvents(session_name,session_onset,session_dur,session_amp);
        session_st.regressor_no_interest = true;
        keys{N_distinctvelocities+1} =  session_name;
        values{N_distinctvelocities+1} = session_st;
        
        % session1_name = 'Session1';
        % session1_onset = onset(1)-70;
        % session1_dur = session_dur;
        % session1_amp = 0.1;
        % session1_st = nirs.design.StimulusEvents(session1_name,session1_onset,session1_dur,session1_amp);
        % session1_st.regressor_no_interest = true;
        % keys{N_distinctvelocities+1} =  session1_name;
        % values{N_distinctvelocities+1} = session1_st;
        % 
        % 
        % 
        % session2_name = 'Session2';
        % session2_onset = onset(end)-70;
        % session2_dur = session_dur;
        % session2_amp = 0.1;
        % session2_st = nirs.design.StimulusEvents(session2_name,session2_onset,session2_dur,session2_amp);
        % session2_st.regressor_no_interest = true;
        % keys{N_distinctvelocities+2} =  session2_name;
        % values{N_distinctvelocities+2} = session2_st;      
    end
    values{ii} = st;
end

if length(unique(keys))~= length(keys) 
    %THIS SHOULD NOT HAPPEN, 
    %unless two velocities are very close to each other and when rounded, 
    %produce the same integer. if it does happen, you can combine the 2 or
    %give them different tags inside this 'if' statement.
    keyboard
    s=dbstack;
    warning([mfilename, ' line ', num2str(s(2).line), ': Two different velocities with the same tag.']) 
end


ALLSTIM = Dictionary(keys,values);
