function [axpos] = subplotpos(nrow,ncol,varargin)
% axpos = SUBPLOTPOS(NROW,NCOL) returns the (normalized) axes positions of nrow by ncol
%     subplot panels . Unlike "subplot", SUBPLOTPOS neither returns an
%     axes handle nor creates it, it only returns an Nby4 array of 
%     positions. The positions can then be used by axes() function. 
%     Compared to subplot, the empy space between subplot panels is 
%     minimal by default and customizable via input parameters 'xgap' and 'ygap'. 
%     Once "axpos" is returned by SUBPLOTPOS, one can modify other properties of an axes such as 'XColor' to make it look even cleaner, e.g. 
%           axes('Position', [axpos(1,:)],'XColor','w','YColor','w','XTick',[],'YTick',[])
%
% axpos = SUBPLOTPOS(NROW,NCOL,Name,Value); 
%     Input 3 to 10 are optional input parameters in Name-Value pairs: 
%     axpositionz = subplotpos(nrow,ncol,'xgap',value,'ygap',value,'lr',value,'bu',value);
%         'xgap' :	the vertical gap length(empty space) in between panels. Unit is normalized
%           to the width of each panel, i.e. xgap=0.1 means 10 percent of each panel's width is empty space. 
%           The gap cannot be less than 0 or greater than 0.5 (50 percent). Default xgap is 0.002;
%         'ygap' :	the horizontal gap length(empty space) in between panels. Unit is normalized
%           to the width of the panel, i.e. ygap=0.4 means 40 percent of each panel's height is empty space. 
%           The gap cannot be less than 0 or greater than 0.5 (50 percent). Default ygap is 0.002;
%         'lr'   :	left-right: x positions (normalized) of far left and far right among all subplot panels, e.g. [0.1 .9].
%           'lr' controls the empty space width of far left and far right of the figure.
%           'lr' default: [0.0800 0.9200] 
%         'bu'   :	bottom-up:  y positions (normalized) of bottom and top of among all subplot panels, e.g. [0.1 .9].         
%           'bu' controls the empty space height of bottom and top of the figure.
%           'bu' default: [0.0790 0.9400]
%
% EXAMPLE 1
%     axpos = subplotpos(3,2)
%     axpos =
%         0.0804    0.6533    0.4192    0.2864
%         0.5004    0.6533    0.4192    0.2864
%         0.0804    0.3663    0.4192    0.2864
%         0.5004    0.3663    0.4192    0.2864
%         0.0804    0.0793    0.4192    0.2864
%         0.5004    0.0793    0.4192    0.2864
%     The output positions can be used by axes() function, by a command as simple as 
%     for i=1:size(axpos,1)
%         ax(i) = axes('Position',axpos(i,:));
%     end
%     See the demo section at the end of the code for another example.
%     
%EXAMPLE 2
%     axpos = subplotpos(3,2,'xgap',0.03,'ygap',0.05,'lr',[0.01 .99],'bu',[0.01 0.99])
%     axpos =
%         0.0174    0.6715    0.4753    0.3103
%         0.5073    0.6715    0.4753    0.3103
%         0.0174    0.3448    0.4753    0.3103
%         0.5073    0.3448    0.4753    0.3103
%         0.0174    0.0182    0.4753    0.3103
%         0.5073    0.0182    0.4753    0.3103
%
%See also AXES, POSITION, SUBPLOT 
%
%mohsen hozan 6/21/2016

%% Parsing and Validating inputs
narginchk(2,10)%nrow and ncol are the two mandatory inputs 
% nargoutchk(0,1)

classes_general     = {'numeric'};
attributes_general  = {'nonempty','nonnegative','real'};

validateattributes(nrow,    classes_general,    [attributes_general,'scalar','integer'])
validateattributes(ncol,    classes_general,    [attributes_general,'scalar','integer'])

npnl = nrow*ncol;   %number of subplot panels
if npnl>1999
    error('You certainly cannot fit %d subplot panels in a single figure, can you? \nSUBPLOTPOS refuses to proceed. retry with smaller nrow/ncol values.',npnl)
end

p = inputParser;
paramnames = {...
    'xgap';         ... %x-gap: the vertical gap lengths in between panels. Unit is normalized to the width of the panel, i.e. xgap=0.1 means 10 percent of each panel's width is empty space. The gap cannot be greater than 50 percent.
  	'ygap';      	...	%y-gap the horizontal gap lengths in between panels. Unit is normalized to the width of the panel, i.e. ygap=0.1 means 10 percent of each panel's height is empty space. The gap cannot be greater than 50 percent.
    'lr';         	... %left-right: x positions (normalized) of far left and far right among all subplot panels, e.g. [0.1 .9].	MATLAB Defaults: [0.1300 0.9050]
    'bu';         	... %bottom-up:  y positions (normalized) of bottom and top of among all subplot panels, e.g. [0.1 .9].         MATLAB Defaults: [0.1100 0.9250]
    };

defaultz = {...
    0.002;              ... %xgap
    0.002;              ... %ygap
    [0.0800 0.9200];    ... %left-right: x positions (normalized) of far left and far right among all subplot panels, e.g. [0.1 .9]
    [0.0790 0.9400];    ... %bottom-up:  y positions (normalized) of bottom and top of among all subplot panels, e.g. [0.1 .9]
    };

validationFns = {...
    @(x) (validateattributes(x,classes_general,[attributes_general,'scalar','<=',0.5]));                    ... %xgap
    @(x) (validateattributes(x,classes_general,[attributes_general,'scalar','<=',0.5]));                    ... %ygap
    @(x) (validateattributes(x,classes_general,[attributes_general,'increasing','numel',2,'<=',1]));     	... %lr
    @(x) (validateattributes(x,classes_general,[attributes_general,'increasing','numel',2,'<=',1]));      	... %bu
    };

for prm = 1:length(paramnames)
    addParameter(p,paramnames{prm},defaultz{prm},validationFns{prm})
end
parse(p,varargin{1:end})

%% 
xgap    = p.Results.xgap;
ygap    = p.Results.ygap;
lr      = p.Results.lr;
bu      = p.Results.bu;

%%
fig_width   = lr(2) - lr(1);            %proportional height of the figure to be taken by all the axes
fig_height  = bu(2) - bu(1);            %proportional width of the figure to be taken by all the axes
pnlwidth   	= (1-xgap)*fig_width/ncol;  %the width of each panel
pnlheight  	= (1-ygap)*fig_height/nrow; %the height of each panel

midpoints_x = linspace(lr(1),lr(2),2*ncol+1);
midpoints_x = midpoints_x(2:2:end);
midpoints_y = linspace(bu(1),bu(2),2*nrow+1);
midpoints_y = midpoints_y(2:2:end);

xoffsets = midpoints_x - pnlwidth/2; %start xpos of each panel (bottom left corner)
xoffsets = repmat(xoffsets.',[nrow,1]);
yoffsets = midpoints_y(end:-1:1) - pnlheight/2; %start ypos of each panel (bottom left corner)
yoffsets = repmat(yoffsets,[ncol,1]);
yoffsets = yoffsets(:); %to be consistent with MATLAB's incosistent manner of addressing subplot panels. 

widthz   = pnlwidth*ones(npnl,1);
heightz  = pnlheight*ones(npnl,1);

%% Output
axpos    = [xoffsets,yoffsets,widthz,heightz]; %output

%% Demo
demo = false; %change this to "true" to see how the function works. change this back to "false" when leaving the function.
if demo
    ax = gobjects(npnl,1); %initialize an array for graphics objects (axes)
    fh = figure(1001); clf;
    set(fh,'name','subploutpos demo','NumberTitle','off')
    N=100;
    x=linspace(0,2,N);
    for i = 1:size(axpos,1)
        ax(i) = axes;
        clr = 0.9*rand(1,3);
        plot(ax(i),x,randn(1,N),'color',clr)
        text(ax(i),0,3,num2str(i,'%04.0f'),'HorizontalAlignment','left','VerticalAlignment','cap','color',clr,'FontUnits','Normalized','FontSize',0.25)
        axis(ax(i),[0 2 -3 3])
        ax(i).Position = axpos(i,:);
        ax(i).XColor = [1 1 1];
        ax(i).YColor = [1 1 1];
        ax(i).Box = 'off';
        ax(i).XTick = [];
        ax(i).YTick = [];
%         drawnow
    end
end

% linkaxes(ax,'x')

