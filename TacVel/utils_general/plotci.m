function [varargout] = plotci(t,x,lowerbound,upperbound,varargin)
% [object_handle] = plotci(t,x,lowerbound,upperbound)
%   Plot x(t) with shaded confidence intervals specified by upperbound and
%   lowerbound
%[object_handle] = plotci(t,x,upperbound,lowerbound,Name,Value)
%   Name-value input pairs are optional.
%   Valid Name-Value pairs include:
%       'color': color of the filled confidence intervals, as specifiec by
%               a colorspec (see "doc colorspec").
%       'alpha': transparency of the filled confidence interval, a positive
%               number between [0 1] inclusive.
%
% EXAMPLE:
%       t=linspace(0,4*pi,1000);
%       x1=sin(t)+.05*randn(size(t));
%       upperbound1 = sin(t)+1;
%       lowerbound1 = 0.75*sin(t)-1;
%       ci1 = plotci(t,x1,upperbound1,lowerbound1,'color',[0 1 0],'alpha',0.3);
%
%       x2=(randn(size(t))./15)+exp(t./max(t))-1.5;
%       upperbound2 = x2+.5;
%       lowerbound2 = x2-.5;
%       hold on
%       ci2 = plotci(t,x2,upperbound2,lowerbound2,'color','cyan','alpha',0.4);
%       hold off
%
% Morgan Siegmann 10/21/2016
%% Parsing and Validating inputs
narginchk(4,8)%t,x,lowerbound and upperbound are mandatory inputs
% optional inputs: color, alpha
classes_general     = {'numeric'};
attributes_general  = {'nonempty','real'};

validateattributes(t,classes_general,[attributes_general,'vector'])
validateattributes(x,classes_general,[attributes_general,'vector'])
validateattributes(upperbound,classes_general,[attributes_general,'vector'])
validateattributes(lowerbound,classes_general,[attributes_general,'vector'])

p = inputParser;
paramnames = {...
    'color';... %color: input in the syntax for specifying color in matlab:RGB triplet,Short name,Long name.
    'alpha' };  %alpha: transparency of the fill color; scalar in the range of [0,1]

defaultz = {...
    'green';... %color
    0.3 };      %alpha

validationFns = {...
    @(x) (validateattributes(x,{'char','numeric'},{}));...
    @(x) (validateattributes(x,{'numeric'},{'scalar','<=',1,'>=',0,}))};

addParameter(p,paramnames{1},defaultz{1},validationFns{1})
addParameter(p,paramnames{2},defaultz{2},validationFns{2})

parse(p,varargin{1:end})

color = p.Results.color;
alpha = p.Results.alpha;

color2 = [0 0 0];
if ~ischar(color)
    color2= 1*color;
else
    % Mohsen hates string colorspecs
end

%%
t=t(:).';
x=x(:).';
lowerbound=lowerbound(:)';
upperbound=upperbound(:)';

X=[t,fliplr(t)];
Y=[lowerbound,fliplr(upperbound)];
%%
ax = gca;
hold on
plot(ax,t,x,'color',color2,'LineWidth',2)
ci = fill(ax,X,Y,color,'FaceAlpha',alpha,'EdgeAlpha',0.2,'EdgeColor',color2);
hold off
if nargout
    varargout{1} = ci;
end