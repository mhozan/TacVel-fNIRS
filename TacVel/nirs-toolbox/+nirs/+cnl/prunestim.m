classdef prunestim < nirs.modules.AbstractModule
    %% PRUNESTIM prunes the Stimulus 
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
    properties
        % % discard stimulus events, e.g. intro, session1, etc.
        % 
        stims2keep = {'intro','v20cmps','v25cmps','v31cmps','v39cmps'};
        % indexstim2keep = zeros(size(stims2keep));
        trial_onset_shift = 0;
        trial_dur = 15;
        unifystims = false; %Default is false. make True if comparing stim-on vs stim-off conditions are desired.
        % % values_stim = raw.stimulus.values;
        % % keys_stim = raw.stimulus.keys;
        % for kk=1:length(stims2keep)
        %     indexstim2keep(kk) = find(contains(raw.stimulus.keys,stims2keep{kk}));
        % end
        % 
        % values_stim = raw.stimulus.values(indexstim2keep);
        % keys_stim = raw.stimulus.keys(indexstim2keep);
        % raw.stimulus = Dictionary(keys_stim,values_stim);
    end

    methods
        function obj = prunestim( prevJob )
            obj.name = 'nirs_cnl_prunestim';

            if nargin > 0
                obj.prevJob = prevJob;
            end
        end

        function data = runThis( obj, data )
            indexstims2keep = zeros(size(obj.stims2keep));
            N_prune = length(indexstims2keep);
            for i = 1:numel(data)                
                if isempty(data(i).stimulus)
                    return
                end                
                
                keys = data(i).stimulus.keys;

                for aa=1:N_prune
                    indexstims2keep(aa) = find(contains(keys,obj.stims2keep{aa},'IgnoreCase',true));
                end

                values_pruned = data(i).stimulus.values(indexstims2keep);
                keys_pruned = keys(indexstims2keep);

                for cc=1:N_prune
                    values_pruned{cc}.onset = values_pruned{cc}.onset + obj.trial_onset_shift;
                    values_pruned{cc}.dur = obj.trial_dur*ones(size(values_pruned{cc}.dur));
                end

                if obj.unifystims && N_prune>1
                    values_unified= values_pruned{1};

                    for jj=2:N_prune
                    values_unified.onset =   [values_unified.onset;values_pruned{jj}.onset];
                    values_unified.dur =   [values_unified.dur;values_pruned{jj}.dur];
                    values_unified.amp =   [values_unified.amp;values_pruned{jj}.amp];
                    end

                    values_unified.name = 'stimON';

                    values_pruned = {values_unified};
                    keys_pruned = {'stimON'};
                end

                stims_pruned = Dictionary(keys_pruned,values_pruned);
                data(i).stimulus = stims_pruned;

            end
        end
    end
end


