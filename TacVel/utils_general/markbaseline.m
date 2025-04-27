function [] = markbaseline(varargin)
% marks the baseline on the x-axis of the all the axes in the given figure(or
% current figure if no input) with a dark vertical-dashed line.
%%%%%%%%Example 1: markbaseline([3600])
%           drwas a marker line at x=3600 on all the axes of the current figure;
%%%%%%%%Example 2: markbaseline([-10 10 100])
%           drwas three marker lines at x=-10 & x=10 & x=1000 on all the axes of the current figure;
%%%%%%%%Example 3: markbaseline(figurehandle,[3600])
%           drwas a marker line at x=3600 on ALL the axes of the figure given by figurehandle;
%%%%%%%%Example 4: markbaseline(axeshandle,[3600])
%           drwas a marker line at x=3600 only on the axes given by axeshandle;
%by Mohsen hozan@mit.edu 11/14/2016
narginchk(1,2)
%
% if nargin == 0
%     fh=gcf;
% else
%     fh=varargin{1};
% end
if isgraphics(varargin{1},'figure')
    fh = varargin{1};
	AllAxes = findobj(fh,'type','axes');
    if nargin~=2
        error('There must be a second input argument containing the x-location(s) of baseline marks.')
    end
elseif isgraphics(varargin{1},'axes')
    fh = gcf;
    AllAxes = varargin{1};
    if nargin~=2
        error('There must be a second input argument containing the x-location(s) of baseline marks.')
    end
elseif isnumeric(varargin{1})
    fh = gcf;
    AllAxes = findobj(fh,'type','axes');
    x = varargin{1};  x=x(:);
    if nargin>1
        error('When the first input argument is numeric, the second input argument is invalid/unnecessary.')
    end
else
	error('Invalid input arguments.')
end
if nargin ==2
    if isnumeric(varargin{2})
        x=varargin{2}; x=x(:);
    else
        error('Second input argument must be numeric.')
    end
end

x = repmat(x,1,2);

% AllAxes = findobj(fh,'type','axes');

for ax=1:length(AllAxes)
    for mark = 1:size(x,1)
    xmark = x(mark,:);
    y = AllAxes(ax).YLim;
    hold(AllAxes(ax),'on')
    line(AllAxes(ax),xmark,y,'color',[.1 .1 .1],'linewidth',1,'linestyle','--')
    end
end
