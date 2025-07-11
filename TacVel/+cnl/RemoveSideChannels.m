classdef RemoveSideChannels < nirs.modules.AbstractModule
    %% Remove Side Channels - removes the following channels primarily covering the temporal lobe due to their movement artifacts caused by the nirx optodes being loose on the sides of the head

    % S16-D18
    % S16-D16
    % S14-D18
    % S1-D2
    % S1-D1
    % S2-D1
    % S4-D4
    % S12-D13
% 

    methods
        function obj = RemoveSideChannels( prevJob )
            obj.name = 'RemoveSideChannels';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            detectors_to_be_excluded =[1,18,4,13];
            % detectors_to_be_excluded =[1,18,4,13,19:26];
            sources_to_be_excluded =[1,16];
            % detectors_to_be_excluded =[1,18,4,13,2,7,11,16,19:26]; %focus on S1 only (primary somatosensory cortex); also exludes short detectors
            % detectors_to_be_excluded =[1,18,4,13,2,7,11,16]; %focus on S1 only (primary somatosensory cortex)
            % sources_to_be_excluded =[1,4,8,12,16];


            for i = 1:numel(data)
                    if(isa(data,'nirs.core.Data'))
                        lst = ismember(data(i).probe.link.source,sources_to_be_excluded)  | ismember(data(i).probe.link.detector,detectors_to_be_excluded);
                        midline_channels = ismember(data(i).probe.link.source,[8,9])  & ismember(data(i).probe.link.detector,[9,10]);
                        % nnz(lst2) %must be 12
                        lst = or(lst,midline_channels);
                        data(i).probe.link(lst,:)=[];
                        data(i).data(:,lst)=[];
                        
                    elseif(isa(data,'nirs.core.ChannelStats'))
                        lst = ismember(data(i).probe.link.source,[1,16])  | ismember(data(i).probe.link.detector,[1,18,4,13]);
                        data(i).probe.link(lst,:)=[];
                        data(i).variables(lst,:)=[];
                        data(i).beta(lst)=[];
                        data(i).covb(lst,:)=[];
                        data(i).covb(:,lst)=[];
                        
                    end
            end
        end
    end
    
end

