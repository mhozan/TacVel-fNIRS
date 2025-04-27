function [varargout]=printALL(varargin)
%PrintALL saves a PNG version of all the currently open figures.
%   Detailed explanation goes here
%   
%   Input Variables:
%
%   Output Variables:
%
%   Example(s):
%   
%   see also 
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2017

narginchk(0,1)
nargoutchk(0,1)

if nargin==1
    validateattributes(varargin{1},{'char'},{'nonempty'})
    pathstr = varargin{1};
else
    pathstr = fullfile(cd,['pngfolder',datestr(now,'-yyyymmdd')]);
end

allFHs = get(0, 'Children');

if isempty(allFHs)
    disp('No figure found to be printed.')
    return
end



for i = 1:length(allFHs)
    fighndl = allFHs(i,1);
    printPNG(fighndl,pathstr)
end
disp(['PNG figures are succesfully saved in : ',pathstr,'\'])
prompt = 'Do you want to open the folder containing the PNG files? N/Y [N]: ';
tic
str = [];
% str = input(prompt,'s');
if any(strcmpi(str,{'Y','Yes','Yes ','Sure','Absolutely','Why not?','Yes please'}))
    winopen(pathstr)
end

if nargout > 0
    varargout{1} = pathstr;    
end