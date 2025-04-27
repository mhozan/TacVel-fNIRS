function raw = loadNIRx2(folder,registerprobe)
% This function loads NIRx data

% <subjid>.wl1  - wavelength #1 data
% <subjid>.wl2  - wavelength #2 data
% <subjid>_config.txt	- config file
% <subjid>.evt	- stimulus events	  (data taken from config file)
% <subjid>_probeInfo.mat - probe file
% <subjid>.tpl -topology file  (data taken from config file)


if(~isfolder(folder))
    folder=fileparts(folder);
end

if(nargin<2)
    registerprobe=true;  % flag to also do 3D registration
end

disp(['Loading : ' folder]);
probe=nirs.io.loadNIRxProbe(folder,registerprobe);


% this is a hack that allowed me to fix a study where the NIRx aquistion had
% been setup wrong and was missing some channels.  NIRx saves all data in
% the wl files (even those not part of the probe).  Creating a SDmask.mat
% file overwrites the info in the NIRx file
if(~isempty(dir(fullfile(folder,'*SDmask.mat'))))
    fi=dir(fullfile(folder,'*SDmask.mat'));
    load(fullfile(folder,fi(1).name));
    disp('loading mask from file');
    info.S_D_Mask=SD_mask;


    [s,d]=find(info.S_D_Mask);


    %% Fix added based on bug from Guilherme A. Zimeo Morais created an issue 2016-08-24
    [s,a]=sort(s); d=d(a);

    link=table(s,d,...
        'VariableNames',{'source','detector'});
    probe.link=link;

end


file = dir(fullfile(folder,'*.hdr'));
if(isempty(file))
    file = dir(fullfile(folder,'*_config.json'));  % new format
    if(isempty(file))
        raw=[];
        return;
    else
        info = parsehdrJSON(fullfile(folder,file(1).name));
    end
else
    info = parsehdr(fullfile(folder,file(1).name));
end


if(isfield(info,'ShortDetectors') && info.Detectors>size(probe.detPos,1))
    % this is an intermediate version of the software when short
    % seperations were first introduced.
    lst=find(ismember(probe.optodes.Type,{'Source'}));
    lst2=find(ismember(probe.optodes.Type,{'Detector'}));
    lst3=find(~ismember(probe.optodes.Type,{'Source','Detector'}));

    DetShort =probe.optodes(lst,:);
    DetShort3D =probe.optodes_registered(lst,:);
    for i=1:length(lst)
        str=['000' num2str(i+info.Detectors-info.ShortDetectors)];
        DetShort.Name{i}=['Detector-' str(end-3:end)];
        DetShort.Type{i}='Detector';
        DetShort.X(i)=DetShort.X(i)+eps(1);
        DetShort3D.Name{i}=['Detector-' str(end-3:end)];
        DetShort3D.Type{i}='Detector';
        DetShort3D.X(i)=DetShort3D.X(i)+1;

    end
    probe.optodes=[probe.optodes(lst2(1:info.Detectors-info.ShortDetectors),:); ...
        DetShort; probe.optodes(lst,:); probe.optodes(lst3,:)];
    probe.optodes_registered=[probe.optodes_registered(lst2(1:info.Detectors-info.ShortDetectors),:); ...
        DetShort3D; probe.optodes_registered(lst,:); probe.optodes_registered(lst3,:)];

    info.S_D_Mask(:,info.Detectors-info.ShortDetectors+1:end)=eye(info.ShortDetectors);


end


%% Now, let's get the data
file = rdir(fullfile(folder,'*.nirs' ));
if(~isempty(file))
    raw=nirs.io.loadDotNirs(file(1).name,true);
    probe.link=raw.probe.link;
    raw.probe=probe;

else
    [s,d]=find(info.S_D_Mask);
    [s,a]=sort(s); d=d(a);

    link=table(s,d,...
        'VariableNames',{'source','detector'});
    if(height(probe.link)<height(link))
        probe.link=link;
    end

    probe.link=[[probe.link; probe.link] ...
        table(reshape(repmat(info.Wavelengths(:)',...
        height(probe.link),1),[],1),'VariableNames',{'type'})];





    % read the standard wl files
    raw = nirs.core.Data();
    raw.probe=probe;

    kk=find(ismember(probe.link.type,probe.link.type(1)));
    lst=find(ismember(info.SDkey(:,2:3),[probe.link.source(kk) probe.link.detector(kk)],'rows'));
    for idx=1:length(info.Wavelengths)
        file = dir(fullfile(folder,['*.wl' num2str(idx)]));

        if(file(1).bytes==0)
            error('zero byte file');
        end

        d = dlmread(fullfile(folder,file(1).name));
        raw.data=[raw.data d(:,lst)];
    end

    raw.time=[0:size(raw.data,1)-1]/info.SamplingRate;
end

raw.probe.fixeddistances=raw.probe.swap_reg.distances;


if(isfield(info,'ShortDetectors') && ~isempty(info.ShortDetectors) && ...
        info.ShortDetectors)

    j=nirs.modules.LabelShortSeperation;
    j.max_distance=10;
    raw=j.run(raw);
    % raw.probe.link(ismember(raw.probe.link.detector,info.ShortDetIndex),:)

end

if(isfield(info,'FileName'))
    raw.description=info.FileName;
else
    raw.description=folder;
end

% Add the demographics info
file = dir(fullfile(folder,'*.inf'));
if(~isempty(file))
    demoinfo = parsehdr(fullfile(folder,file(1).name));
else
    file = dir(fullfile(folder,'*_description.json'));
    if(~isempty(file))
        demoinfo = nirs.io.loadjson(fullfile(folder,file(1).name));
    end
end

%% add some of the header (hdr) information to the last field of the demographics field
% this is not the perfect place to put this information, but i prefer it to
% having to modify the nirs.core.data function to allow for properly named
% additional fields.
fldnames=fieldnames(info);
for ff=1:length(fldnames)
    % demoinfo = setfield(demoinfo,fldnames{ff},getfield(info,fldnames{ff}));
    demoinfo.(fldnames{ff})=info.(fldnames{ff});
end
% demoinfo.Header_info = info;
%% add the excel spreadsheet info to demographics field
[xlsdata,xlsinfo] = readExcelCsvFile( oneupdirectory(folder,0));
if isempty(xlsdata) %check one upper folder
    [xlsdata,xlsinfo] = readExcelCsvFile(oneupdirectory(folder,1));
end
%% neatify and integrate the excel info into  raw.demographics
if isempty(xlsdata)
    %Do nothing

else
    demoinfo.Gender = upper(demoinfo.Gender);

    warning off
    xlsdata3 = rows2vars(xlsdata);
    warning on

    AgeYears = xlsdata3.Age_years_Months_;
    demoinfo.AgePrecise = str2double(AgeYears{1});
    if numel(AgeYears)>1 %month is entered
        demoinfo.AgePrecise = demoinfo.AgePrecise+ AgeYears{2}/12;
        xlsdata3(2,:)=[]; %delete the second column entirely
    end
    demoinfo.DateOfExperiment = xlsdata3.DateOfExperiment{:};
    demoinfo.Initials = upper(xlsdata3.Initials{:});
    demoinfo.HandMeas_cm_1to2 = str2double(xlsdata3.HandMeasurement_1_2_);
    demoinfo.HandMeas_cm_2to3 = str2double(xlsdata3.HandMeasurement_2_3_);
    demoinfo.HandMeas_cm_3to4 = str2double(xlsdata3.HandMeasurement_3_4_);
    demoinfo.HandMeas_cm_4to5 = str2double(xlsdata3.HandMeasurement_4_5_);
    demoinfo.HandMeas_cm_1to5 = demoinfo.HandMeas_cm_1to2+demoinfo.HandMeas_cm_2to3+demoinfo.HandMeas_cm_3to4+demoinfo.HandMeas_cm_4to5;
    demoinfo.HandDominant = upper(xlsdata3.Handedness_L_R_{:});
    demoinfo.HairColor = lower(xlsdata3.HairColor{:}(1:5));
    demoinfo.HairType = lower(xlsdata3.HairType{:});
    demoinfo.HairLength = lower(xlsdata3.HairLength{:});
    demoinfo.HeadCircumference_mm = 10*str2double(xlsdata3.HeadCircumference_cm_); %nirs.modules.Resize_Heads takes this by default
    demoinfo.HeadSize = str2double(xlsdata3.HeadCircumference_cm_);
    demoinfo.Capsize = str2double(xlsdata3.CapSize);
    demoinfo.NotesonChannels = xlsdata3.fnirsYellowOrRedChannels{:};
    demoinfo.GalileoSeqBeginning = xlsdata3.GalileoSequences{:};
    demoinfo.NotesOnExperiment=[(xlsdata(end,1).Row{:}),'. ',(xlsdata(end-1,1).Row{:}),'. ',(xlsdata(end-2,1).Row{:}),'.'];
    % contains(xlsdata.Row,'Row')
end

%% make some binary(additional) categorical variables for experimentation with the mixed effects model
demoinfo.isMale = any(strcmpi(demoinfo.Gender,{'Male','M'}));
demoinfo.isCapTight = demoinfo.HeadSize-demoinfo.Capsize>0; %A loose cap is prone to movement artifacts, a tight cap is uncomfortable and affects ECBF
demoinfo.isHairShort = strcmpi(demoinfo.HairLength,'Short'); %Longer hair is bulkier, adversely affects optode-skin contact
demoinfo.isHairLight = any(strcmpi(demoinfo.HairColor,{'Blonde','Blond','White'}));
demoinfo.isHairStraight = any(strcmpi(demoinfo.HairType,{'straight'})); %curly or way hair adversely affects optode-skin contact
demoinfo.isRightHanded = any(strcmpi(demoinfo.HandDominant,{'Right','R'}));
demoinfo.isHandBig = demoinfo.HandMeas_cm_1to5>20.5; % median hand measurements 1_5 across 26 subjects is 20.5

%% folder info placeholder
demoinfo.isFirstSession = contains(folder,'Session1','IgnoreCase',true);
demoinfo.ID = demoinfo.Name; %placeholder; will be filled in nirs.cnl.loadDirectory2 in line 82;
if demoinfo.isFirstSession
    demoinfo.Session = '1'; %placeholder; will be filled in nirs.cnl.loadDirectory2 in line 82
elseif contains(folder,'Session2','IgnoreCase',true)
    demoinfo.Session = '2'; %placeholder; will be filled in nirs.cnl.loadDirectory2 in line 82
end

%% make the variables categorical for statistical analysis
demoinfo.Gender = categorical({demoinfo.Gender}); %make it categorical
demoinfo.HairColor = categorical({demoinfo.HairColor}); %make it categorical
demoinfo.HairType = categorical({demoinfo.HairType}); %make it categorical
demoinfo.HairLength = categorical({demoinfo.HairLength}); %make it categorical

%% remove some demographics info
demo=Dictionary;
chandisstr = num2str(demoinfo.ChanDis,'%2.0f,');
demoinfo.ChanDis = chandisstr(1:end-1);
fieldstoremove = {'APD',...usually empty; causes error in nirs.modules.mixedmodel
    'DateOfExperiment',...Date of experiment is the same as Date
    'Gains',... 2D matrix; causes an error in nirs.modules.MixedEffects
    'Events',... 2D matrix; causes an error in nirs.modules.MixedEffects
    'S_D_Mask',... 2D matrix; causes an error in nirs.modules.MixedEffects
    'Wavelength1',... 2D matrix; causes an error in nirs.modules.MixedEffects
    'Wavelength2',... 2D matrix; causes an error in nirs.modules.MixedEffects
    'TrigIns',...
    'TrigOuts',...
    'AnIns',...
    'Threshold',...
    };

demoinfo = removefields(demoinfo,fieldstoremove);
flds=fields(demoinfo);

for idx=1:length(flds)
    demo(flds{idx})=demoinfo.(flds{idx});
end

raw.demographics=demo;


% There is a slight error in the NIRx files for hyperscanning that I am
% going to exploit to add this info.  For hyperscanning files, the
% info.S_D_Mask covers only 1 subject but the info.SD_Key field is the full
% (both subjects) length.
if(isfield(info,'ChanDis'))
    if(length(info.ChanDis)==height(probe.link))
        raw.demographics('hyperscan')=info.FileName;
    end
end


% Now add stimulus info
if(isfield(info,'Events') && ~isempty(info.Events))
    stimName = unique(info.Events(:,2));
    stimulus=Dictionary();
    for idx=1:length(stimName)
        s = nirs.design.StimulusEvents();
        s.name=['channel_' num2str(stimName(idx))];
        s.onset=info.Events(find(info.Events(:,2)==stimName(idx)),1);
        s.dur=ones(size(s.onset));
        s.amp=ones(size(s.onset));
        stimulus(s.name)=s;
    end
    raw.stimulus=stimulus;
end


if(exist(fullfile(folder,'stimulus.mat')))
    load(fullfile(folder,'stimulus.mat'))
    raw.auxillary('stim')=raw.stimulus;
    raw.stimulus=stimulus;
    disp('loading stim-events from file');
end


if(~isempty(rdir(fullfile(folder,['*_nirsInfo.mat']))))
    f=rdir(fullfile(folder,['*_nirsInfo.mat']));
    info=load(f(1).name);
    if(isfield(info.nirsInfo,'eventInfo'))
        for i=1:length(info.nirsInfo.eventInfo.events.names)
            st=nirs.design.StimulusEvents;
            st.name=info.nirsInfo.eventInfo.events.names{i};
            st.onset=info.nirsInfo.eventInfo.events.onsets{i};
            st.dur=info.nirsInfo.eventInfo.events.durations{i};
            st.amp=ones(size(st.onset));
            raw.stimulus(st.name)=st;
        end
    end
end






%% Add the ADI data (.adicht) as auxillary
loadADICHT = true;
trimtimevectors = true; %flag to trim the excess end of the ADIdata and the beginnig of the fNIRS data
% headerinfoindex= find(contains(raw.demographics.keys,'Header_info'));
startDATEindex = find(matches(raw.demographics.keys,'Date'));
startTIMEindex = find(matches(raw.demographics.keys,'Time'));

nirx_recordstart_str =  ...
    ([raw.demographics.values{startDATEindex}(6:end),' ',...
    raw.demographics.values{startTIMEindex}]); %use for synchronization of ADI and NIRScout

% nirx_recordstart_str =  ...
%     ([raw.demographics.values{headerinfoindex}.Date(6:end),' ',...
%     raw.demographics.values{headerinfoindex}.Time]); %use for synchronization of ADI and NIRScout
stimulus_delay = 0;
if ~isempty(raw.stimulus.values)
    stimulus_delay = raw.stimulus.values{1}.onset(1);
end
ADIdata = [];
if loadADICHT
    ADIdata = nirs.cnl.loadADICHT( oneupdirectory(folder,0),nirx_recordstart_str,stimulus_delay);
    if isempty(ADIdata)
        ADIdata = nirs.cnl.loadADICHT( oneupdirectory(folder,1),nirx_recordstart_str,stimulus_delay);
    end
end

if ~isempty(ADIdata)
    if ~isempty(raw.auxillary.values)
        warning([info.FileName,...
            ': Overwriting the already existing auxillary data by the following ADI dataset:  ',...
            ADIdata.values{1}.extra.ADICHTfilepath])
    end
    if trimtimevectors
        tmin = max(raw.time(1),ADIdata.values{1}.time(1));
        N_todiscard_atthebeginning = nnz(raw.time<tmin);
        raw.time(1:N_todiscard_atthebeginning) = [];
        raw.data(1:N_todiscard_atthebeginning,:) = [];
    end


    raw.auxillary = ADIdata;
end



%% replace the stimulus with renamed eventized triggers using nirs.cnl.eventizeGtriggers2
experiment2021protocol = raw.stimulus.count>0;
cleanupviamask = true;%experiment2021protocol;
resampleaux = false;

% If the deepest directroy folder names are either Session1 or Session2, that means I have duplicated the data into two folder so that I can zero out the second half of the Session1 folder and zero out the first half of the Session2 folder.
% Sessionize = strcmpi(folder(end-8:end-2),'Session') ;
Sessionize = contains(folder,'Session','IgnoreCase',false) ;
replacewith = 0;
if Sessionize
    stimvec= raw.stimulus.values{1}.onset;
    stimdiff = diff(stimvec);
    [~,intermission_start_index] = max(stimdiff);
    intermission = stimvec(intermission_start_index:intermission_start_index+1);
    % halfpoint = floor(size(raw.data,1)/2); %This does not work for 20211025

    % intermission = find(pwr);
    extrapadding = 100; %seconds
    halfpoint_1 = intermission(1)+extrapadding;
    [~,halfpoint_1_index] = min(abs(raw.time-halfpoint_1));
    halfpoint_2 = intermission(2)-extrapadding;
    [~,halfpoint_2_index] = min(abs(raw.time-halfpoint_2));

    % figure; subplot(2,1,1), raw.draw;
    % sgtitle(raw.demographics.values{1})
    % subplot(2,1,2), plot(stimdiff)
    if contains(folder,'Session1')
        % raw.data(halfpoint:end,:) = replacewith;
        raw.data(halfpoint_1_index:end,:) = replacewith;
        % secondhalf = find(raw.stimulus.values{1}.onset>halftime);
        % raw.stimulus.values{1}.amp(secondhalf) = zeros(size(secondhalf));
        % subplot(2,1,2), raw.draw
        % title('Session1')
    elseif contains(folder,'Session2')
        % raw.data(1:halfpoint,:) = replacewith;
        raw.data(1:halfpoint_2_index,:) = replacewith;
        % subplot(2,1,2), raw.draw
        % title('Session2')
    else
        keyboard
        %do nothin
    end
end



if experiment2021protocol && any(contains(raw.auxillary.keys,'Galileo'))
    % raw = raw_all(18);
    % raw2 = raw_all(18);
    onset_shift = stimulus_delay;%raw.stimulus.values{1}.onset(1);%+raw.auxillary.values{1}.record_start_delay_s;
    galileoindex = find(contains(raw.auxillary.keys,'Galileo')); %usually number 3
    triggerAnalog = raw.auxillary.values{galileoindex}.data;
    % pwr = movmean(sum((triggerAnalog),2).^2,1000);
    % figure; subplot(3,1,1), raw.draw; hold on, plot(pwr)
    % subplot(3,1,2), plot(raw.time,pwr)
    triggerDigitized = triggerAnalog>(max(triggerAnalog)/2);
    triggerDigitized_onsets = (find(triggerDigitized));
    triggerDigitized_onsetdiffs = diff(triggerDigitized_onsets);
    [~,triggerDigitized_diffmax_onset]= max(triggerDigitized_onsetdiffs);
    intermission = triggerDigitized_onsets(triggerDigitized_diffmax_onset:triggerDigitized_diffmax_onset+1);
    extrapadding2 = 10000; %samples
    halfpointTrig_1 = intermission(1)+extrapadding2;
    halfpointTrig_2 = intermission(2)-extrapadding2;
    % figure, clf,
    % subplot(5,1,1), plot(triggerAnalog),
    % sgtitle(raw.demographics.values{1})
    if Sessionize
        % halfpoint2 = floor(length(triggerAnalog)/2); %This does not work for 20211025

        if contains(folder,'Session1')
            % triggerAnalog(halfpoint2:end)=replacewith;
            triggerAnalog(halfpointTrig_1:end)=replacewith;
        elseif contains(folder,'Session2')
            % triggerAnalog(1:halfpoint2)=replacewith;
            triggerAnalog(1:halfpointTrig_2)=replacewith;
        else
            %should not happen
            keyboard
        end
    end
    % subplot(5,1,2), plot(triggerDigitized),
    % subplot(5,1,3), plot(triggerDigitized_onsets)
    % subplot(5,1,4), plot(triggerDigitized_onsetdiffs)
    % subplot(5,1,5), plot(triggerAnalog)

    minimum_known_velocity = 20;  % For the 2021 dataset, 20cmps is the lowest velocity.
    Fs_aux=raw.auxillary.values{galileoindex}.Fs; %1000 Hz is typically the sampling frequency of ADI
    trial_onset_shift = 0;  % %stim_on duration is 15 seconds within each 40s trial; this forwards the onset of each trial for analysis by 15s.
    trial_dur = 20; %stim_off duration is 25 seconds within each 40s trial; if the trial_onset_shift is greater than 0, this should be less than 40, so that the sum is never greater than 40 to avoid overlapping with the subsequent trial.
    if demoinfo.isFirstSession && strcmp(demoinfo.Name,'20221130ADM')
    %there is a single solitary lone isolated weird trigger on this
    %subjects which will be ignored by the following command. 
    triggerAnalog(1622360:1622395)=.09;
    disp(mfilename)
    disp([demoinfo.Name , ' Session ', demoinfo.Session, ': triggerAnalog issue is fixed. Line 464.'])
    end

    stimulus = nirs.cnl.eventizeGtriggers2(triggerAnalog,onset_shift,minimum_known_velocity,Fs_aux,trial_onset_shift,trial_dur);

    for ii =1:raw.stimulus.count %preserve original stimulus data as a regressor of no interest.
        stimulus_original.values{1,ii} = raw.stimulus.values{ii};
        stimulus_original.values{1,ii}.regressor_no_interest = true;
        stimulus_original.keys{1,ii} = raw.stimulus.keys{ii};
    end
    keys = horzcat(stimulus.keys,stimulus_original.keys);
    values = horzcat(stimulus.values,stimulus_original.values);
    stimulus_combined = Dictionary(keys,values);
    if Sessionize %do not preserve the Channel1 stimulus
        stimulus_combined = Dictionary(stimulus.keys,stimulus.values);
    end
    raw.stimulus = stimulus_combined; %replace the stimulus with enhanced tagged events



    % figure(101), clf, subplot(2,1,1), raw.draw; subplot(2,1,2), raw2.draw;
    if cleanupviamask && any(contains(raw.stimulus.keys,'Session','IgnoreCase',true))
        % keyboard
        Sessionindex = find(contains(raw.stimulus.keys,'Session'));

        s_onset = raw.stimulus.values{Sessionindex}.onset;
        s_finish = s_onset+raw.stimulus.values{Sessionindex}.dur;

        tvec1 = raw.time;
        mask_raw = ones(size(tvec1));
        mask_raw(tvec1<s_onset) = replacewith;
        mask_raw(tvec1>s_finish) = replacewith;
        % mask_raw(tvec1<s2onset & tvec1>s_finish) = replacewith;
        % % % Session1index = find(contains(raw.stimulus.keys,'Session1'));
        % % % Session2index = find(contains(raw.stimulus.keys,'Session2'));
        % % % 
        % % % s1onset = raw.stimulus.values{Session1index}.onset;
        % % % s1finish = s1onset+raw.stimulus.values{Session1index}.dur;
        % % % s2onset = raw.stimulus.values{Session2index}.onset;
        % % % s2finish = s2onset+raw.stimulus.values{Session2index}.dur;
        % % % 
        % % % tvec1 = raw.time;
        % % % mask_raw = ones(size(tvec1));
        % % % mask_raw(tvec1<s1onset) = replacewith;
        % % % mask_raw(tvec1>s2finish) = replacewith;
        % % % mask_raw(tvec1<s2onset & tvec1>s1finish) = replacewith;
        % figure; plot(tvec1,2.5*mask_raw)

        data = raw.data;
        raw.data = data.*mask_raw;
        % figur([raw.demographics.subject,'_cleaned']) , clf, raw.draw;

        values_aux = ADIdata.values;
        keys_aux = ADIdata.keys;

        for jj=1:ADIdata.count
            tvec_aux = values_aux{jj}.time;
            mask_aux = ones(size(tvec_aux));
            mask_aux(tvec_aux<s_onset) = replacewith;
            mask_aux(tvec_aux>s_finish) = replacewith;
            % % % mask_aux(tvec_aux>s2finish) = replacewith;
            % % % mask_aux(tvec_aux<s2onset & tvec_aux>s_finish) = replacewith;
            % hold on; plot(tvec_aux,2.5*mask_aux)
            data_aux = values_aux{jj}.data;
            values_aux{jj}.data = data_aux.*mask_aux;
        end
    end

end

end

%% parsehdrJSON
function info =parsehdrJSON(file)
% This sub-routine parses the NIRx header info

info=nirs.io.loadjson(file);

info.S_D_Mask=[];
for i=1:length(info.channel_mask)
    for j=1:length(info.channel_mask{i})
        info.S_D_Mask(i,j)=str2num(info.channel_mask{i}(j));
    end
end

info.det_Mask=[];
for i=1:length(info.det_split)
    for j=1:length(info.det_split{i})
        info.det_Mask(i,j)=str2num(info.det_split{i}(j));
    end
end
lst=find(sum(info.det_Mask)==8);
if(~isempty(lst))
    info.ShortDetIndex=find(info.det_Mask(:,lst));
    info.ShortDetectors=length(info.ShortDetIndex);
else
    info.ShortDetIndex=[];
end
info.Detectors=length(sum(info.det_Mask));

info.Wavelengths=[760 850];
info.Sources=length(info.drv_amplitudes);

info.SDkey=[1:info.Sources*info.Detectors;repmat(1:info.Sources,1,info.Detectors); repmat(1:info.Detectors,1,info.Sources)]';

end

%% parsehdr
function info =parsehdr(file)
% This sub-routine parses the NIRx header info

info=struct;
fid=fopen(file,'r');

while(1)
    line=fgetl(fid);
    if(~isstr(line))
        break
    end
    if(contains(line,'='))
        if(~contains(line,'="#'))
            fld = line(1:strfind(line,'=')-1);
            val = line(strfind(line,'=')+1:end);
            if(~contains(val,'"'))
                val=str2num(val);
            else
                val(strfind(val,'"'))=[];
            end
        else
            % matrix- read until next #
            fld = line(1:strfind(line,'=')-1);
            cnt=1; val=[];
            while(1)
                line=fgetl(fid);
                if(~contains(line(1),'#'))
                    val(cnt,:)=str2num(line);
                    cnt=cnt+1;
                else
                    break;
                end
            end
        end
        fld(strfind(fld,'-'))='_';
        fld(strfind(fld,' '))='_';
        info=setfield(info,fld,val);
    elseif(contains(line,'['))
        %header skip
    end

end
fclose(fid);

if(isfield(info,'Wavelengths'))

    % fix a few strings that should be numeric
    info.Wavelengths=str2num(info.Wavelengths);
    info.Mod_Amp=str2num(info.Mod_Amp);
    info.Threshold=str2num(info.Threshold);
    info.ChanDis=str2num(info.ChanDis);

    % Fix the SD key
    if(strcmp(info.S_D_Key(end),',')); info.S_D_Key(end)=[]; end
    keys=strsplit(info.S_D_Key,',');
    for idx=1:length(keys)
        a=strsplit(keys{idx},{'-',':'});
        info.SDkey(idx,1)=str2num(a{3});  % channel index
        info.SDkey(idx,2)=str2num(a{1});  % source
        info.SDkey(idx,3)=str2num(a{2});  % detector
    end
end

try
    % added to fix bug issue #57.  11/1/2018
    if isfield(info, 'ShortDetIndex') && isfield(info,'ShortBundles') && info.ShortBundles > 0
        info.ShortDetectors = length(str2num(info.ShortDetIndex));
    end
end

end