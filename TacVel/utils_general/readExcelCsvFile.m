function [varargout] = readExcelCsvFile(folderPath)
%readExcelCsvFile takes a folder path as input, which may contain one 
% Excel or CSV file.
% It reads the file and returns its data and details in a structure.
%   
%   Input Variables:
%
%   Output Variables:
%
%   Example(s):
%   
% Communication Neuroscience Laboratories,
% Center for Brain, Biology & Behavior
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2023
% 
%   see also 
% 
nargoutchk(1,2)
% Get a list of all files and folders in the input folder
contents = dir(folderPath);

% Remove directories from the list of contents
contents = contents(~[contents.isdir]);

% Filter the list to only include files with the desired extensions
fileList = {contents(~[contents.isdir] & ...
    (contains({contents.name},'.xls') | ...
    contains({contents.name},'.xlsx') | ...
    contains({contents.name},'.csv'))).name};

fileList(contains(fileList,'~$'))=[]; %ignore temporary locked files by microsoft office

% Check if the folder contains only one Excel or CSV file
if numel(fileList) ~= 1
    % error('The input folder must contain exactly one Excel or CSV file.');
    disp('The input folder must contain exactly one Excel or CSV file.');
    varargout{1}=[];
    if nargout==2
        varargout{2}=[];
    end
    return
end

% Get the file name and extension
[~, fileName, fileExtension] = fileparts(fileList{1});

% Read the file
if strcmpi(fileExtension, '.csv')
    data = readtable(fullfile(folderPath, fileList{1}), 'Delimiter', ',');
else
    data = readtable(fullfile(folderPath, fileList{1}),'ReadRowNames',true,'NumHeaderLines',0);
end
% data = rmmissing(data); %remove missing




% Get file details
fileDetails.fileName = fileName;
fileDetails.fileExtension = fileExtension;
fileDetails.filepath = folderPath;
fileDetails.numRows = size(data, 1);
fileDetails.numCols = size(data, 2);
fileDetails.headers = data.Properties.VariableNames;

varargout{1}=data;
if nargout > 1
    varargout{2}=fileDetails;
end 

end

