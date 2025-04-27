function [s_synced]=syncstructs(s,newFs,fname2,fname1,varargin)
% Use this function to compensate for the different sampling frequencies in
% the DEX project data
%
% Morgan Siegmann
% December, 02, 2016
%% Input Arguments
% inputs:
% s: structure
% field names to be synchronized (Y) - in timeseries data it is your
%       signal - want to make it accept more than one of these
% field names to be synchronized (X) - in timeseries data it is the
%       time vector
% new sampling frequency
% Field names must be strings that follow matlab variable naming rules
%validateattributes
validateattributes(s,{'struct'},{})
validateattributes(newFs,{'numeric'},{'integer'})
validateattributes(fname1,{'char'},{})
validateattributes(fname2,{'char'},{})

nfield = numel(varargin)+1;
nfieldname = [{fname1},varargin];

%%
% disp([mfilename, ': Adding synchronized fields to the structure and saving the results...'])
% savedresults = '\\mickey\MIT\Morgan\DEX\DEX Delta\deltapower.mat';
% load(savedresults);
% z=who('-file',savedresults)

SampFreq=cell2mat({s.Freqs});

[uniqfreqs,ind]= unique(SampFreq);
tvec1 = s(ind(1)).(fname2);
tvec2 = s(ind(2)).(fname2);

% if there are no additional input arguments do this once
%   else: do the same thing again but replace fname1 with fname3,4,5 ect.
for i = 1:nfield
    
    x = {s.(nfieldname{i})};
    indcs250 = SampFreq==uniqfreqs(1);
    indcs254 = SampFreq==uniqfreqs(2);
    
    xf250 =  cell2mat(x(indcs250).');
    xf254 =  cell2mat(x(indcs254).');
    if isvector(xf250)
        xf250 =  cell2mat(x(indcs250)).';
        xf254 =  cell2mat(x(indcs254)).';
    end
    
    [tvec_common,data1_synced,data2_synced] = synchronize2(tvec1,xf250,tvec2,xf254,newFs);
    % downsamp_factor=newFs/10;
    % tvec_common = downsample(tvec_common,downsamp_factor);
    % data1_synced = downsample(data1_synced.',downsamp_factor).';
    % data2_synced = downsample(data2_synced.',downsamp_factor).';
    
    data_synced = zeros(size([data1_synced;data2_synced]));
    data_synced(indcs250,:) = data1_synced;
    data_synced(indcs254,:) = data2_synced;
    
    %% Remove the 'time vector' field ***%make this flexible for multiple fields
    s=rmfield(s,nfieldname{i});
    
    %% Repopulate the Structure ***%make this flexible for multiple fields
    for n = 1:length(s)
        s(n).(nfieldname{i}) = data_synced(n,:);
        %     s(i).tvec_common = tvec_common.';
        %     disp(i)
    end
end
s=rmfield(s,{fname2});% only need to do fname2 once (not included in extra fields feature)
s(1).(fname2) = tvec_common;
%% Keep this
s_synced=s;
%% Save Results
% newfilename=[savedresults(1:end-4),'_synced','.mat'];
% save(newfilename,'s');
%
% disp('New file is saved at:')
% disp(newfilename)
