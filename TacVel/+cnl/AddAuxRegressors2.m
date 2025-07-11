classdef AddAuxRegressors2 < nirs.modules.AbstractModule
    %% AddAuxRegressors2 - Adds auxillary data as regressors to the GLM model
    % edited by hozan@huskers.unl.edu
    
    properties
        normalize;  % normalize the regressors
        orth;
        label;
    end
    
    methods
        function obj = AddAuxRegressors2( prevJob )
            obj.name = mfilename;
            obj.normalize = true;
            obj.label={'mic';'EMG';'ADI';'accel';'Resp';'Eyeblink';'NIRmidline';'iblink';'hand'};
            %             obj.label={'axis'};
            %             obj.label={'EMG'};
            obj.orth = false;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                
                auxlabels=data(i).auxillary.keys;
                %                 idx=find(ismember(lower(auxlabels),lower(obj.label)));
                idx = contains(auxlabels,obj.label,'IgnoreCase',true);
                if(ismember('all',lower(obj.label)))
                    idx=logical(1:data(i).auxillary.count);
                end
                for j=1:length(obj.label)
                    if(~isempty(strfind(obj.label{j},'*')))
                        ll=strsplit(obj.label{j},'*');
                        for k=1:length(auxlabels)
                            found=true;
                            for l=1:length(ll)
                                if(~isempty(ll{l}))
                                    found=found & ~isempty(strfind(auxlabels{k},ll{l}));
                                end
                            end
                            if(found)
                                idx=[idx k];
                            end
                        end
                    end
                end
                
                if(isempty(idx))
                    continue;
                end
                aux = data(i).auxillary.values(idx);
                if(obj.orth)
                    dd=[];                    
                    for j=1:length(aux)
                        for k=1:size(aux{j}.data,2)
                            %                             dd=[dd interp1(aux{j}.time,aux{j}.data(:,k),data(i).time)];
                            dd=[dd interp1(aux{j}.time,aux{j}.data(:,k),data(i).time,'linear','extrap')]; %default linear interpolation produces NaN
                        end
                    end
                    if(obj.normalize)
                        dd=dd-ones(size(dd,1),1)*mean(dd,1);
                        dd=dd./(ones(size(dd,1),1)*std(dd,[],1));
                    end
                    dd=orth(dd); %toolbox\matlab\matfun\orth.m
                    for j=1:size(dd,2)
                        st=nirs.design.StimulusVector;
                        st.regressor_no_interest=true;
                        st.name=['AuxPCA' num2str(j)];
                        st.time=data(i).time;
                        st.vector=dd(:,j);
                        data(i).stimulus(st.name)=st;
                    end
                else
                    for j=1:length(aux)
                        dd=[];
                        
                        for k=1:size(aux{j}.data,2)
                            %                             dd=[dd interp1(aux{j}.time,aux{j}.data(:,k),data(i).time)];
                            dd=[dd interp1(aux{j}.time,aux{j}.data(:,k),data(i).time,'linear','extrap')]; %default linear interpolation produces NaN
                        end
                        
                        if(obj.normalize)
                            if any(isnan(dd))
                                keyboard
                            end
                            dd=dd-ones(size(dd,1),1)*mean(dd,1);
                            if std(dd,[],1) %avoid creating an INF matrix by dividing by 0
                                dd=dd./(ones(size(dd,1),1)*std(dd,[],1));
                            end
                        end
                        % dd(isnan(dd)) = 0; %take care of NANs
                        if any(isnan(dd))
                            keyboard
                        end
                        dd=orth(dd);
                        for k=1:size(dd,2)
                            st=nirs.design.StimulusVector;
                            st.regressor_no_interest=true; %default:true 
                            %                             st.name=[data(i).auxillary.keys{idx(j)} num2str(k)];
                            st.name=aux{j}.description;
                            st.time=data(i).time;
                            st.vector=dd(:,k);
                            data(i).stimulus(st.name)=st;
                        end
                    end
                end
                
                
                
                
            end
        end
    end
    
end

