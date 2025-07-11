function [varargout] = adi2auxnirs(ADICHTfilepath,varargin)
%ADI2AUXNIRS creates the "auxillary" field formatted for the Huppert lab 
%nirs toolbox. The input is the full path to the *.adicht file.
%
%   Example(s):
%   auxillary = adi2auxnirs('Z:\folder\sample.adicht');
%
% Written by:
% Mohsen Hozan hozan@huskers.unl.edu
% Communication Neuroscience Laboratories
% Center for Brain, Biology, and Behavior
% University of Nebraska-Lincoln
% November 2020
%
% see also
% nirs, Dictionary, nirs.design.StimulusEvents, nirs.cnl.ADIdata
nargoutchk(0,1)
narginchk(0,3)

standardize_zscore = false; %centered to have mean 0 and scaled to have standard deviation 1
[~,adifilename,~]=fileparts(ADICHTfilepath);

%make sure ADI_matlab_sdk for reading .adicht files is on the matlab path
if ~contains(path,'JimHokanson-adinstruments_sdk_matlab','IgnoreCase',true)
    addpath('Z:\students\mhozan2\Matlab Libraries\JimHokanson-adinstruments_sdk_matlab-1061f17')
end

%make sure Huppert nirs toolbox is on the matlab path
if ~contains(path,'nirs-toolbox','IgnoreCase',true)
    addpath(genpath('Z:\students\mhozan2\Matlab Libraries\nirs-toolbox-master-20250416\'))    
end


file_header = adi.sdk.openFile(ADICHTfilepath);
adidatafileinfo = adi.readFile(ADICHTfilepath);
Nch = adidatafileinfo.n_channels; %number of channels (typically 9)
% Nrec =  adidatafileinfo.n_records %number of records (typically 1); change to 1 if errors occur
Nrec =  1;

recorddurations = zeros(adidatafileinfo.n_records,1);
for rec = 1:adidatafileinfo.n_records
recorddurations(rec) = adidatafileinfo.records(rec).duration;
end
[~,indmaxdurrec] = max(recorddurations); %find the index of the record with maximum duration


% chnames = adidatafileinfo.channel_names;
isDOUBLEprecision = true; %'single' precision if false, 'double' if true.
keys = adidatafileinfo.channel_names; %channel names
for ki = 1:length(keys) %makes sure channel names are unique, and valid variable names, otherwise Dictionary() throws an error.
    keys{ki} = matlab.lang.makeValidName([num2str(ki,'CH%02d_'),keys{ki}]);  
end

values = cell(Nrec,Nch);
% for  ii =1:1
%     name = char(velocitylabels(ii));
%     keys{i_rec,i_ch} = name;
%     onset = trainonsets(stimvelocityorder==ii);
%     dur = 10*ones(size(onset)); %seconds
%     amp = 1*ones(size(onset));
%     st = nirs.design.StimulusEvents(name,onset,dur,amp);
%     
%     values{ii} = st;
% end

%%
   

for i_rec=indmaxdurrec% 1:Nrec
    for i_ch=1:Nch
        ch_info =  adidatafileinfo.channel_specs(i_ch);
        
        if ch_info.n_samples(i_rec) == 0
            values{1,i_ch} = [];
            continue
        end
     	description = keys{1,i_ch};

        data  = adi.sdk.getChannelData(file_header,...
            i_rec,...record id
            ch_info.id,...channel id
            1,...first sample
            ch_info.n_samples(i_rec),...%number of samples
            isDOUBLEprecision);
        data = data(:); %columnize for nirs toolbox compatibility
%         isnormal = kstest(data);
%         isnormal = lillietest(data,'Alpha',0.01);
%         isnormal = false;
%         str = [];
%         if isnormal
%             str = sprintf([mfilename, ' :     Normal :',adifilename,' : ', description,' is Normal, i.e. comes from a standard normal distribution.\n']);
%         else
%             str = sprintf([mfilename, ' : Not Normal :',adifilename,' : ', description,'  is NOT Normal, i.e. does not come from a standard normal distribution.\n']);
%         end
%         fprintf(str);
        
        if standardize_zscore
            data= zscore(data);
        end
        
        link = table();  %if multiple timeseries are within the same 
        %channel, nirs toolbox links the names and indices with variable.
        time = (0:ch_info.n_samples(i_rec)-1).'/ch_info.fs(i_rec);
%         Fs = ch_info.fs; %dependent variable
        
        %the following field are extra info and not interprable by the nirs
        %toolbox's default methods.
        extra.ADICHTfilepath = ADICHTfilepath;
        extra.downsample_amount = ch_info.downsample_amount(i_rec);
        extra.trigger_minus_rec_start_seconds = ... %this variable is pretty crucial in synchronizing the NIRX and ADI signals.
            adidatafileinfo.records(i_rec).trigger_minus_rec_start;
        extra.trigger_minus_rec_start_samples = ...
            floor(extra.trigger_minus_rec_start_seconds*ch_info.fs(i_rec));
        extra.data_starts = ch_info.data_starts(i_rec);
        extra.record_starts = ch_info.record_starts(i_rec);
        extra.comments = adidatafileinfo.records(i_rec).comments;
        extra.units = ch_info.units{i_rec};
        extra.data_start_str = adidatafileinfo.records(i_rec).data_start_str;

        %%calculate ADI delay and add it to the adidata
        record_start_delay_sec = 0; %default value
        discardthedelay = false;
        if nargin == 2
            nirx_recordstart  = datevec(varargin{1}{1});
            adi_recordstart = datevec(adidatafileinfo.channel_specs(1).record_starts(i_rec));
            record_start_delay_vec = adi_recordstart - nirx_recordstart; %delay datetime vector
            record_start_delay_sec = record_start_delay_vec*datetimeconversion2seconds;

            stimulus_delay = varargin{1}{2};
            if stimulus_delay>0
            discardthedelay = true;
            end
            
        end


        if record_start_delay_sec < 0
            disp(record_start_delay_sec)
% keyboard
        end
        if datevec(adidatafileinfo.channel_specs(1).record_starts(i_rec)) ~= datevec(adidatafileinfo.channel_specs(1).data_starts(i_rec))
         
        end
        values{1,i_ch} = nirs.cnl.ADIdata(data,time,link,description,extra,record_start_delay_sec,stimulus_delay, discardthedelay);
        

        
    end
end
%remove channels with no recorded samples
containssamples = ~cellfun('isempty',values);
values = values(containssamples);
keys = keys(containssamples);

% discardindices = contains(keys,'galileo','IgnoreCase',true) | ...
%     contains(keys,'Hi_Resp','IgnoreCase',true) |...
%     contains(keys,'AbLo','IgnoreCase',true) |...
%     contains(keys,'xaccl','IgnoreCase',true) | ...
%     contains(keys,'yaccl','IgnoreCase',true);
% 
% indices_to_be_kept = ~discardindices;find(~contains(keys,'galileo','IgnoreCase',true));
% values = values(indices_to_be_kept);
% keys = keys(indices_to_be_kept);

auxillary = Dictionary(keys,values);

if nargout ==1
    varargout{1} = auxillary;
end

