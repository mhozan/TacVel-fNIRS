function [ hrate_bpm ,varargout] = heartbeatestimate_freqbased(ECG_vec,varargin)
% HEARTBEAT_ESTIMATE estimates the heart rate in beats per minute and
% returns  it in a timeseries of the same size as the input.
%
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2021
% See also
% nirs.muscle.estimate_heartrate, hannbinning

%% 
fs = 1000;
if nargin>=2
    fs = varargin{1};
end
%%
windur_s = 10; %seconds
prcntoverlap = .9;
heartbeatrange = [45 115]/60; %heartbeats per minute (divided by 60 to be in Hz)

% nfft = 2^(3+nextpow2(winLen));
nfft = 2^16;

%%
winLen = windur_s*fs;
noverlap = prcntoverlap*winLen;
winStep = winLen - noverlap;


%Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
ECG = hannbinning(ECG_vec,winLen,prcntoverlap);
ECG = fft(ECG,nfft);
ECG = abs(ECG/nfft);
ECG = ECG(1:nfft/2+1,:);
ECG(2:end-1,:) = 2*ECG(2:end-1,:);
f_vec = fs*(0:(nfft/2))/nfft;
%%
L = size(ECG,2);
tvec2 = (0:1:L-1).';
tvec2 = tvec2*winStep/fs ;

timeshift = 0;
tvec2 = tvec2 + windur_s/2 +timeshift;

%%


[~,fmin_index] = min(abs(f_vec-heartbeatrange(1)));
[~,fmax_index] = min(abs(f_vec-heartbeatrange(2)));
f_vec = f_vec(fmin_index:fmax_index);
ECG = ECG(fmin_index:fmax_index,:);
bpm_vec = 60*f_vec;


[~,heartrate_freqindex_Hz] = max(ECG);
hrate_bpm = 60*f_vec(heartrate_freqindex_Hz);
%%
% beatsperminute = hrate_bpm;
% tvec2 =linspace(0,);
figur(mfilename),clf,
ax=axes;
imagesc(ax,tvec2,bpm_vec,ECG) 
ax.Position= subplotpos(1,1,'bu',[.4,.97]);
ax.YDir = 'normal';
% ax.YLim = ylimz*60;
% ax.XLim = xlimz;
% ax.CLim = climz;
title('Spectrogram-Extracted Heartbeat')
xlabel('t (S)')
ylabel('heart beats per minute')
% cl = colorbar;
% cl.Label.String = '|P1(f)|';
% ('Label','Hi');

hold on 
plot(ax,tvec2,hrate_bpm,'r','LineWidth',2)

% ax2= axes;
% plot(ax2,tvec2,ECG_vec)
% % ax2.XLim = xlimz;
% ax2.Position= subplotpos(1,1,'bu',[.05,.4]);
% linkaxes([ax,ax2],'x')

%%
if nargout>1
    varargout{1} = [];
end
