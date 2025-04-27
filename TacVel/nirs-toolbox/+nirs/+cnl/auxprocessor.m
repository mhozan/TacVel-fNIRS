function [nirsdata] = auxprocessor(obj,varargin)
%AUXPROCESSOR is a custom script that selects certain auxillary channels
%and discards the rest. It also splits and processes the EMG channels to
%two distinct high frequency and low frequency (moving average) channels.
%It labels the regressor_no_interest flag of each aux signal according to
%default settings. feel free to modify these settings inside the file.
%   Example(s):
%    [nirsdata] = auxprocessor(nirsdata,varargin);
%
% Written by:
% Mohsen Hozan hozan@huskers.unl.edu
% Communication Neuroscience Laboratories
% Center for Brain, Biology, and Behavior
% University of Nebraska-Lincoln
% Apr 2021
%
% see also
% nirs, Dictionary, nirs.design.StimulusEvents, nirs.cnl.adi2auxnirs


%%
aux = obj.auxillary;
N_aux=aux.count;
values = aux.values;
keys = aux.keys;
doignorecase = true;

%% tag which labels to keep and which to ignore/discard
auxlabelstokeep = {'mic';'EMG';'axis';'ADI'};

for ii = 1:N_aux  %% select certain aux channels and discard the rest
    if ~contains(keys{ii},auxlabelstokeep,'IgnoreCase',doignorecase)
        values{ii}.keepit = false;
    end
end


%%

%% resample
newFs = 1000;
% obj.Fs 
for jj=1:N_aux
    if values{jj}.Fs == newFs
        continue %skip
    end
        
    % resample data
    d = values{jj}.data;
    t = values{jj}.time;
    
    % de-mean the data to avoid edge effects
    mu = nanmean(d);
    d = bsxfun(@minus,d,mu);
    
    [d,new_t]=resample(d,t,newFs);
    % restore original mean
    d = bsxfun(@plus,d,mu);
    
%     values{jj}.Fs = newFs;
   	values{jj}.data = d;
    values{jj}.time = new_t;
	values{jj}.resampled = true;
end
%% accelerometer low pass and detrend
for xx=1:N_aux
    if ~contains(keys{xx},'axis','IgnoreCase',doignorecase)
        continue %skip
    end
    
    %low pass filter at 20 Hz
    f_lowpass = 20; %hz
    values{xx}.data = lowpass_accelerometer(values{xx}.data,values{xx}.Fs,f_lowpass);   
    values{xx}.description = [values{xx}.description,'_lpf',num2str(f_lowpass)];
    keys{xx} = [keys{xx},'_lpf',num2str(f_lowpass)];
    
    %linear piece-wise detrend every 30 seconds
%     values{xx}.data = detrend_lin_pw(values{xx}.data,30*values{xx}.Fs); 
%     values{xx}.description = [values{xx}.description,'_lindetrend'];
%     keys{xx} = [keys{xx},'_lindetrend',num2str(f_lowpass)];  
%     
%     figure(3874), clf,
%     ax(1)=subplot(2,1,1);plot(values{xx}.time,values{xx}.data)
%     ax(2)=subplot(2,1,2);plot(aux.values{xx}.time,aux.values{xx}.data)
%     linkaxes(ax,'x')

end

%% ADI pulse and heartrate
windur_s = 6; %seconds
prcntoverlap = .9;


for yy=1:N_aux
    if ~contains(keys{yy},'ADIpulse','IgnoreCase',doignorecase)
        continue %skip
    end
%     winLen = windur_s*values{yy}.Fs;
%     winStep = winLen - noverlap;
%     noverlap = prcntoverlap*winLen;
%     % nfft = 2^(3+nextpow2(winLen));
%     nfft = 2^16;    

    ECG_vec = values{yy}.data;

    ECG = hannbinning(ECG_vec,winLen,prcntoverlap);
ECG = fft(ECG,nfft);

L = size(ECG,2);
tvec2 = (0:1:L-1).';
tvec2 = tvec2*winStep/fs ;
tvec2 = tvec2 + windur_s/2 +tvec(1);
%Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.

ECG = abs(ECG/nfft);
ECG = ECG(1:nfft/2+1,:);
ECG(2:end-1,:) = 2*ECG(2:end-1,:);
f = fs*(0:(nfft/2))/nfft;
end


% ekg_analysis20210326
%nirs.muscle.estimate_heartrate
%% split EMGs to high freq and low freq (mvgavg)
f_hipass = 50;
count= 0; %counter for adding processed EMG signals

for zz=1:N_aux
    if ~contains(keys{zz},'EMG','IgnoreCase',doignorecase)
        continue %skip
    end
    disp(keys{zz})
    
    count = count+2;
    newkey_mvgavg = [keys{zz},'_movingavg200ms'];
    newkey_hpf = [keys{zz},'_hipass50hz'];
    
    
    
    %make 2 copies
    values{N_aux+count-1} = values{zz}; 
    keys{N_aux+count-1} = newkey_mvgavg; 
    values{N_aux+count}   = values{zz};
    keys{N_aux+count} = newkey_hpf; 
    
    %change the keepit flag of the original copy to false, i.e. discard it later
	values{zz}.keepit = false; 

    %mvg avg
    data_aux_lp = movingavg_abs_emg(values{zz}.data,values{zz}.Fs);%-nanmean(data_aux);
    values{N_aux+count-1}.data = data_aux_lp;
    values{N_aux+count-1}.description = newkey_mvgavg;
    
    %hi pass
    data_aux_hp = highpass_emg(values{zz}.data,values{zz}.Fs,f_hipass);
 	values{N_aux+count}.data = data_aux_hp;
    values{N_aux+count}.description = newkey_hpf;
    
    
%     figure(7863), clf, 
%     ax(1)=subplot(3,1,1); plot(values{zz}.time,values{zz}.data), title('raw')
%     ax(2)=subplot(3,1,2); plot(values{zz}.time,data_aux_hp), title('highpass')
%     ax(3)=subplot(3,1,3); plot(values{zz}.time,data_aux_lp), title('lowpass')
%     linkaxes(ax,'x')

end



%%
tvec_aux = eval(varnamez{vv}).auxillary.values{aux_count}.time;
data_aux = eval(varnamez{vv}).auxillary.values{aux_count}.data;
data_aux = zscore(data_aux);
if eval(varnamez{vv}).auxillary.values{aux_count}.extra.downsample_amount == 1
    tvec_aux = downsample(eval(varnamez{vv}).auxillary.values{aux_count}.time,40);
    data_aux = downsample(eval(varnamez{vv}).auxillary.values{aux_count}.data,40);
    %                     keyboard
end




fs = 1000;

if contains(key,{'EMG'},'IgnoreCase',true)  %highpass
    f_hipass = 50;
    data_aux2 = highpass(data_aux,f_hipass,fs,'ImpulseResponse','fir');
    data_aux1 = movingavg_abs_emg(data_aux,fs);%-nanmean(data_aux);
    data_aux = [data_aux2,data_aux1];
    
end


if contains(key,{'axis'},'IgnoreCase',true)
    f_lowpass = 40;
    data_aux = lowpass(data_aux,f_lowpass,fs)-nanmean(data_aux);
end
if contains(key,{'ADI'},'IgnoreCase',true)
    f_lowpass = 80;
    data_aux = lowpass(data_aux,f_lowpass,fs)-nanmean(data_aux);
end



%%


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

nirsdata = Dictionary(keys,values);
