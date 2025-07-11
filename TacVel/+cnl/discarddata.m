function nirsobject_shortened = discarddata(nirsobject,T_s,varargin)
%DISCARDDATA deletes any data before a given time instant (in seconds)
%   RAW_d = DISCARDDATA(RAW,T_s) returnes a shortened version of RAW
%   with discarded data and stimulus values before T_s.
%
%   RAW_d = DISCARDDATA(RAW,T_s,'isdicardbefore',true) also returnes a shortened version
%   of RAW with discarded data and stimulus values before T_s. The last
%   input paramter is an optional binary with true for discarding before
%   T_s and false for discarding after T_s.
%
%
%   RAW_d = DISCARDDATA(RAW,T_s,'isdicardbefore',false) returnes a shortened version
%   of RAW with discarded data and stimulus values after T_s.
%
%   Input Variables:
%       nirsobject : nirs class data with 'stimulus', 'time', and 'data'
%       attributes. Its size can be 1x1 or 1xn
%       T_s : time instant in seconds, where data has to be discarded
%     	T_s can be a scalar or 1xn vector matching the size of nirsobject
%   Optional Input Pair (logical):
%       isdicardbefore: true (default) for discard before T_s or false
%       for discard after T_s. isdicardbefore must be a 1x1 or 1xn
%       logical array matching the size of nirsobject
%   Output Variables:
%       nirsobject_shortened : a shortened version of RAW with appropriate deletion
%
%   Example(s):
%   RAW_discardedbefore_20s  = nirs.cnl.discarddata(RAW,20,'before');
%   RAW_discardedafter_3000s = nirs.cnl.discarddata(RAW,3000,'after');
%
% Written by:
% Mohsen Hozan hozan@huskers.unl.edu
% Communication Neuroscience Laboratories
% Center for Brain, Biology, and Behavior
% University of Nebraska-Lincoln
% November 2020
%
% see also
% nirs
%% Parsing and Validating inputs
narginchk(2,4)  %nirsobject and T_s are the two mandatory inputs
nargoutchk(0,1) % if no output arg, the input nirs data will be modified.

classes_general     = 'numeric';
attributes_general  = 'nonempty';

validateattributes(nirsobject,    {'nirs.core.Data'},...
    {attributes_general})

validateattributes(T_s,     {classes_general},    ...
    {attributes_general,'nonnegative','vector'})
T_s = T_s(:); %columnize T_s
p = inputParser;
paramnames = {...
    'isdiscardbeforeTs';	... %binary input. if false, data after T_s will be discarded.
    };
defaultz = {...
    true;              ... %default discarding happens to data before T_s
    };
validationFns = {...
    @(x) (validateattributes(x,{'numeric','logical'},{attributes_general,'vector','nonnegative','integer'}));
    };

for prm = 1:length(paramnames)
    addParameter(p,paramnames{prm},defaultz{prm},validationFns{prm})
end
parse(p,varargin{1:end})


isbefore    = logical(p.Results.isdiscardbeforeTs(:)); %columnized

%%
L = size(nirsobject,1);
%repmat the single value 2nd and 3rd inputs L times to match the size of nirsobject.

if and(size(T_s,1) ~= L     ,   size(T_s,1) > 1)
    error('Size mismatch! T_s must either be a scalar or a vector as long as the nirs data.')
else
    T_s = repmat(T_s,L,1);  %
end

if and(size(isbefore,1) ~= L    ,   size(isbefore,1) > 1)
    error('Size mismatch! ''isdicardbefore'' must be either 1x1 or the length of the nirs data.')
else
    isbefore = repmat(isbefore,L,1);  %
end

%%
discardindex = zeros(size(T_s));
discardnegativetimes = true; %stimulus and auxillary time vectors can become negative after being shifted. change this flag to discard the negative times.
% discardindex_aux = zeros(size(T_s),nirsobject); %auxillay signals have independent time vectors with different sampling frequencies
for ii=1:L
    [~,discardindex(ii)] = min(abs(nirsobject(ii).time - T_s(ii)));
%     discardindex(ii) = discardindex(ii);
    if isbefore(ii) %discarding before T_s
        
        nirsobject(ii).time(1:discardindex(ii)-1) = [];
        nirsobject(ii).data(1:discardindex(ii)-1,:) = [];
        
    else %discarding after T_s
        nirsobject(ii).time(discardindex(ii)+1:end) = [];
        nirsobject(ii).data(discardindex(ii)+1:end,:) = [];
    end
    
    
    timeshift = nirsobject(ii).time(1);
    nirsobject(ii).time = nirsobject(ii).time - timeshift; %start time at zero
    
    % shift stimulus; discard negative times
    stimname = unique(nirs.getStimNames(nirsobject(ii)));
    if(~iscellstr(stimname)); stimname={stimname}; end
    timeshift_stimz = timeshift*ones(size(stimname));
    for jj = 1:length(stimname)
        st = nirsobject(ii).stimulus(stimname{jj});
        st.onset = st.onset-timeshift_stimz(jj);
        
        if discardnegativetimes
            n_negativeonsets = length(find(st.onset<0));
            st.onset(1:n_negativeonsets)=[];
            st.dur(1:n_negativeonsets)=[];
            st.amp(1:n_negativeonsets)=[];
        end


        nirsobject(ii).stimulus(stimname{jj})=st;
    end
    
    % shift auxillary time vectors
    values = nirsobject(ii).auxillary.values;
    for aa=1:nirsobject(ii).auxillary.count
        tvec_aux = values{aa}.time;
        data_aux = values{aa}.data;
        tvec_aux = tvec_aux-timeshift;
        if discardnegativetimes
            n_negativetimez = length(find(tvec_aux<0));
            tvec_aux(1:n_negativetimez)=[];
            data_aux(1:n_negativetimez)=[];
        end
        values{aa}.time = tvec_aux;
        values{aa}.data = data_aux;
    end
    auxillary = Dictionary(nirsobject(ii).auxillary.keys,values);
    
    nirsobject(ii).auxillary = auxillary;
    %     nirsobject(ii)=nirs.design.shift_stimulus_onset(nirsobject(ii),stimname,shift);
    %     nirsobject(ii).stimulus.va
    % 	nirsobject(i)=nirs.design.shift_stimulus_onset(nirsobject(i),stimname,shift);
end
nirsobject_shortened = nirsobject;



