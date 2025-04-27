function [oneupdir] = oneupdirectory(path,howmanyup)
%ONEUPDIRECTORY returns the path that is 'howmanyup' folder above the 'path'.
%mohsen hozan
%20201111
if nargin==1
    howmanyup=1; %how many folders up in the path?
elseif nargin ==2
    %howmanyup is given
end

validateattributes(path,    {'char'},    {'nonempty','scalartext'})
validateattributes(howmanyup,    {'numeric'},    {'nonnegative','real','scalar','integer'})
[fPath, ~, ~] = fileparts(path);

splitdir = regexp(fPath,filesep,'split');

oneupdir = fullfile(splitdir{1:end-howmanyup});

if isempty(oneupdir)
    error('Too many folders up!')
end
if ~contains(oneupdir,':')  %in case a non-mapped network path is given, which starts with \\ in Windows or // in Mac/Linux
    oneupdir = [filesep,filesep,oneupdir];
end


end

