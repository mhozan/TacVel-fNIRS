classdef DiscardStims < nirs.modules.AbstractModule
%% DiscardStims - Removes specified stim conditions from design.
% 
% Options: 
%     listOfStims - cell array of stim names

    properties
        listOfStims = {}; % cell array of stim names
    end
    
    methods
        function obj = DiscardStims( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Discard Stim Conditions';
        end
        
        function data = runThis( obj, data )
            
           % get current stim names
            stimNames = nirs.getStimNames(data);
            disc=[];
            % loop through and remove stims
            for i = 1:numel(data)
                if(isa(data(i),'nirs.core.Data') | isa(data(i),'eeg.core.Data') | isa(data(i),'nirs.core.GenericData'))
                    data(i).stimulus = data(i).stimulus.remove( obj.listOfStims );
                elseif(isa(data(i),'nirs.core.sFCStats'))
                    cond=data(i).conditions;
                    lst=find(~ismember(cond,obj.listOfStims));
                    data(i).R=data(i).R(:,:,lst);
                    
                    if(length(data(i).dfe)>1)
                        data(i).dfe=data(i).dfe(lst);
                    end
                    data(i).conditions={data(i).conditions{lst}};
                    if(~isempty(data(i).ZstdErr))
                        lst=find(ismember(cond,obj.listOfStims));
                        data(i).ZstdErr(:,:,lst,:)=[];
                        data(i).ZstdErr(:,:,:,lst)=[];
                        
                    end
                else
                    cond=data(i).conditions;
                    cond={cond{~ismember(cond,obj.listOfStims)}};
                    if(~isempty(cond))
                         data(i)=data(i).ttest(cond);
                    else
                        disc=[disc; i];
                    end
                end
            end
            data(disc)=[];
            
            % compare old and new list of stims
            if length(stimNames) == length(nirs.getStimNames(data))
                warning( 'No stims discarded.  Did you provide the correct stim names?' )
            end
        end
    end
    
end

