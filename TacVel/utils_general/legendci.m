function varargout=legendci(varargin)
%legendci('str1','str2',...) add legend to the Line plots drawn by
%'plotci' and ignores adding legend to patches.
% see also
% PLOTCI, LEGEND
% Mohsen Hozan 10/24/2016
linez = findall(gca,'Type','Line'); %finds the LINE objects only.
uistack(linez,'bottom') %brings the LINE objects to the bottom of the stack, for the legend to 

lgnd=legend(varargin);

if nargin~=size(linez,1)
    warning('legendci might be adding legends to the objects other than ''LINE'', which it is not meant to.')
end
if nargout==1
    varargout{1}=lgnd;
end
% set(lgnd,'Interpreter', 'none')