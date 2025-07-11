function ADIdata = loadADICHT( adifiledir , varargin)
%wrapper for handling both directory and file inputs for the adi2auxnirs function, which loads .adicht files
fileExt  = '.adicht';


if isfolder(adifiledir) %input is the folder containing the adicht file
%     files = dir(fullfile(adifiledir,'**',['*',fileExt]));  %** wildcard looks within subfolders
  	files = dir(fullfile(adifiledir,['*',fileExt]));  %** wildcard looks within subfolders
elseif isfile(adifiledir) %input is the full path to the .adicht file.
    files = dir(adifiledir);
else
    error('unknown input file path directory')
end

nFiles = length(files);
% assert(nFiles == 1)
ADIdata = [];

if nFiles>0
ADIdata = nirs.cnl.adi2auxnirs(fullfile(files.folder,files.name),varargin);
end




