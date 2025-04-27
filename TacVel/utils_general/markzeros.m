function [] = markzeros(varargin)
% mark the zeros on the x-axis of the all the axes in the given figure(or
% current figure if no input) with a red vertical line.
%by Mohsen hozan@mit.edu 9/2/2016
narginchk(0,1)

if nargin == 0
    fh=gcf;
else
    fh=varargin{1};
end

AllAxes = findobj(fh,'type','axes');

for ax=1:length(AllAxes)
    x = [0 0];
    y = AllAxes(ax).YLim;
    line(AllAxes(ax),x,y,'linewidth',2','color','r')
end
