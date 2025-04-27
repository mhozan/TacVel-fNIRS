function [fr]=firingrate(timevectorinsecond)

classes     = {'double'};
attributes  = {'nonempty','real'};

validateattributes(timevectorinsecond,    classes,    attributes);
timediff = diff(timevectorinsecond);
fr = 1./timediff;
if any(isinf(fr))
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
    histogram(fr,edges);
    ax=gca;
    set(ax,'XScale','log')
    xlim([edgemin,edgemax])
    % set(ax,'XTick',xtickz,'XTickLabel',xticklabelz);
    xlabel(ax,'Firing Rate (Hz)')
    ylabel(ax,'Spike Count per bin')
    
end


