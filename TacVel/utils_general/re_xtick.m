function []=re_xtick(varargin)
narginchk(0,1)
% set(gcf,'WindowButtonDownFcn','disp(''axis callback'')')
if nargin == 0    
hAllAxes = findobj(gcf,'type','axes');
else
hAllAxes= varargin{1};
end

if ~isempty(hAllAxes)
    for ax=1:length(hAllAxes)        
        ca=hAllAxes(ax);
        xlimz=ca.XLim;
        Xstart=ceil(xlimz(1));
        Xstep = 1;
        Xend = ceil(xlimz(2));
        if range(xlimz)>10*60
            Xstep = 60;
            Xstart = Xstart-mod(Xstart,60);
        end            
        xtickz = Xstart:Xstep:Xend;

        while length(xtickz)>30;  xtickz(2:2:end)=[]; end
        xticklabelz = hms_string(xtickz);
        if length(xticklabelz)>15; xticklabelz(2:2:end,:)=' '; end
        
        set(ca,'XTick',xtickz)        
        if ~isempty(ca.XTickLabel)
            set(ca,'XTickLabel',xticklabelz)
        end
    end
    drawnow
    
else
    %no axes was detected.
end
