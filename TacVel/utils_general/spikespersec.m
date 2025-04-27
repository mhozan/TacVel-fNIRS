function [spkpersec]=spikespersec(spiketimestamps_us,continous_tvec_us)
% SPIKESPERSEC gets the timestamps of the spikes(action potential) in
% microseconds in addition to the continous time vector of the object (in
% us) and returns a vector the same size as the continous time vector with
% the firing rate per one second. 
% mohsen 9/2016
classes     = {'double'};
attributes  = {'nonempty','real'};
validateattributes(spiketimestamps_us,      classes,    attributes);
validateattributes(continous_tvec_us,       classes,    attributes);
Fs_continous = 1e6./mode(diff(continous_tvec_us));
onesec=zeros(Fs_continous,1);
% onesec=zeros(1000,1);

continous_tvec_us = continous_tvec_us(:);
spiketimestamps_us = spiketimestamps_us(:);
spkpersec = zeros(size(continous_tvec_us));
%%
tol= 2000; %tolerance in microsecond (the continous timevec is )
tol = tol./max(abs([spiketimestamps_us;continous_tvec_us]));
[x,y]=ismembertol(spiketimestamps_us,continous_tvec_us,tol);
sum(x)
%%
timediff = diff(spiketimestamps_us);
spkpersec = 1./timediff;
if any(isinf(spkpersec))
    warndlg('Inf in firing rate.')
end
% fr(isinf(fr))= 0;

%%
if true
    edgemax = 10000;
    edgemin = 0.01;
    edges = edgemax;
    Nbinsperpower10tick = 128;
    while edges(1)>edgemin-1e-5
        edges = [edges(1)./10.^(1/Nbinsperpower10tick);edges];
    end
    xtickz = edges;
    xticklabelz = cell(size(edges));
    decimaltickz = find(rem(edges,0.01)<1e-6);
    xticklabelz(decimaltickz)=cellstr(num2str(edges(decimaltickz)));
    figure(321);clf;
    histogram(spkpersec,edges);
    ax=gca;
    set(ax,'XScale','log')
    xlim([edgemin,edgemax])
    % set(ax,'XTick',xtickz,'XTickLabel',xticklabelz);
    xlabel(ax,'Firing Rate (Hz)')
    ylabel(ax,'Spike Count per bin')
    
end


