function [spkpersec]=spikespersec2(spiketimestamps_us,Fs_continous)
% SPIKESPERSEC gets the timestamps of the spikes(action potential) in
% microseconds and returns a time series sampled at Fs with the firing rate per one second. 
% mohsen 9/2016


classes     = {'double'};
attributes  = {'nonempty','real'};
validateattributes(spiketimestamps_us,      classes,    attributes);

continous_tvec_us = spiketimestamps_us(1):1:spiketimestamps_us(end); %out of memory
% Fs_continous = 1e6./mode(diff(continous_tvec_us));
onesec=zeros(Fs_continous,1);

spiketimestamps_us = spiketimestamps_us(:);
spkpersec = zeros(size(spiketimestamps_us));
%%
% tol= 2000; %tolerance in microsecond 
tol = 1e6/Fs_continous/2; %(the continous timevec is usually sampled at around 250 Hz);
tol = eps+tol./max(abs([spiketimestamps_us;continous_tvec_us]));
[x,y]=ismembertol(spiketimestamps_us,continous_tvec_us,tol);
if sum(x)<length(x)
    warning([mfilename,'.m : some spike timestamps are not attributed to the continous time vector.'])
end

time_matches = continous_tvec_us()



%%
