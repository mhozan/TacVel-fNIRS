classdef pruneaux2 < nirs.modules.AbstractModule
    %% PRUNEAUX prunes the aux channels
    %   This module processes various auxiliary signals recorded during 
    %   Near-Infrared Spectroscopy (NIRS) experiments, filters unnecessary 
    %   signals, applies preprocessing (e.g., filtering, thresholding, 
    %   normalization), and retains only relevant auxiliary data.
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
        aux2keep = {'ADIpulsetransducer', 'mic', 'Eyeblink_rate', 'EMG_face', ...
            'hand', 'XYZ_acceleration', 'Respiration_rate', ...
            'CombinedAuxMotion', 'NIRmidline'};
        resampleaux = true;
        trimtimevectors = true;
        clip_threshold = 99;
    end

    methods
        function obj = pruneaux2(prevJob)
            obj.name = 'nirs_cnl_pruneaux';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end

        function data = runThis(obj, data)
            for i = 1:numel(data)
                if isempty(data(i).auxillary)
                    return
                end

                keys_aux = data(i).auxillary.keys;
                values_aux = data(i).auxillary.values;

                % Process different auxiliary signals
                values_aux = obj.processAccelerometer(keys_aux, values_aux);
                values_aux = obj.processRespiration(keys_aux, values_aux, i);
                values_aux = obj.processEyeEMG(keys_aux, values_aux);
                values_aux = obj.processCombinedMotion(keys_aux, values_aux);
                values_aux = obj.processNIRMidline(keys_aux, values_aux, data(i));

                % Keep only selected auxiliary signals
                values_aux = obj.keepSelectedAuxiliary(keys_aux, values_aux);

                % Resample auxiliary data if required
                if obj.resampleaux
                    values_aux = obj.resampleAuxiliary(data(i), values_aux);
                end

                % Store processed auxiliary data
                aux_pruned = Dictionary(keys_aux, values_aux);
                data(i).auxillary = aux_pruned;
            end
        end
    end

    methods (Access = private)
        function values_aux = processAccelerometer(obj, keys_aux, values_aux)
            acclindices = find(contains(keys_aux, 'accl', 'IgnoreCase', true));
            if isempty(acclindices), return; end

            combinedaccelerometer = values_aux(acclindices);
            values_aux{end+1} = combinedaccelerometer{1}; % Duplicate X-accl
            values_aux{end}.description = 'XYZ_acceleration';

            % Apply thresholding and normalization
            for j = 1:3
                data = combinedaccelerometer{j}.data;
                lowerprctile = prctile(data, 90);
                upper_threshold = lowerprctile * 10;

                data = min(max(data, lowerprctile), upper_threshold);
                data = (data - lowerprctile) / max(data - lowerprctile);
                data(isnan(data)) = 1;

                combinedaccelerometer{j}.data = data;
            end

            % Compute vector magnitude (Euclidean norm)
            xyzaccl = sqrt(sum(cellfun(@(x) x.data .^ 2, combinedaccelerometer, 'UniformOutput', false), 2));
            xyzaccl = xyzaccl - median(xyzaccl);

            % Apply RMS envelope filter
            xyz_acc_env = obj.applyFilter(xyzaccl, combinedaccelerometer{1}.Fs, 'rms');

            values_aux{end}.data = xyz_acc_env;
        end

        function values_aux = processRespiration(obj, keys_aux, values_aux, i)
            respindices = find(contains(keys_aux, {'resp', 'cage', 'abdomen'}, 'IgnoreCase', true));
            if isempty(respindices), return; end

            RESP_aux = values_aux(respindices);
            RESP_aux_combined = cell2mat(cellfun(@(x) x.data, RESP_aux, 'UniformOutput', false));

            respnormfactor = prctile(RESP_aux_combined, 90);
            if all(respnormfactor)
                RESP_aux_combined = RESP_aux_combined ./ respnormfactor;
            else
                RESP_aux_combined = zeros(size(RESP_aux_combined));
                disp(['i = ', num2str(i, '%0.2d'), ' Respiration rate data is not good. Set to zero.']);
            end

            RESP_median = median(RESP_aux_combined, 2) / 60;
            values_aux{end+1} = RESP_aux{1}; % Duplicate
            values_aux{end}.description = 'Respiration_rate';
            values_aux{end}.data = RESP_median;
        end

        function values_aux = processEyeEMG(obj, keys_aux, values_aux)
            eyeEMGindex = find(contains(keys_aux, 'eye', 'IgnoreCase', true));
            if isempty(eyeEMGindex), return; end

            eyeEMGdata = values_aux{eyeEMGindex}.data;
            eyeblink = eyeEMGdata(:, 1);
            facialEMG = eyeEMGdata(:, 2) .^ 2;

            % Apply high-pass filtering
            Fs = values_aux{eyeEMGindex}.Fs;
            facialEMG = obj.applyFilter(facialEMG, Fs, 'highpass', 5);
            facialEMG = facialEMG .^ 2;

            values_aux{end+1} = values_aux{eyeEMGindex};
            values_aux{end}.description = 'Eyeblink_rate';
            values_aux{end}.data = eyeblink;

            values_aux{end+1} = values_aux{eyeEMGindex};
            values_aux{end}.description = 'EMG_face';
            values_aux{end}.data = facialEMG;
        end

        function values_aux = processCombinedMotion(obj, keys_aux, values_aux)
            emghand = obj.getAuxData(keys_aux, values_aux, 'hand');
            micdata = obj.getAuxData(keys_aux, values_aux, 'mic');

            if isempty(emghand) || isempty(micdata), return; end

            % Apply filtering
            micdata = obj.applyFilter(micdata, values_aux{1}.Fs, 'highpass', 10);

            % Compute combined motion signal
            combinedmotion = emghand + micdata;
            combinedmotion = combinedmotion .^ 2;
            combinedmotion = obj.applyFilter(combinedmotion, values_aux{1}.Fs, 'highpass', 10);
            combinedmotion = abs(combinedmotion);
            combinedmotion = min(combinedmotion, prctile(combinedmotion, 90));

            values_aux{end+1} = values_aux{1};
            values_aux{end}.description = 'CombinedAuxMotion';
            values_aux{end}.data = combinedmotion;
        end

        function values_aux = processNIRMidline(obj, keys_aux, values_aux, data)
            lst = ismember(data.probe.link.source, [8, 9]) & ismember(data.probe.link.detector, [9, 10]);
            if isempty(lst) || ~nnz(lst), return; end

            midline = median(data.data(:, lst), 2);
            midline = log(midline ./ mean(midline, 'omitnan'));
            midline(midline < 0) = 0;
            midline = midline - mean(midline, 'omitnan');

            values_aux{end+1} = values_aux{1};
            values_aux{end}.description = 'NIRmidline';
            values_aux{end}.data = midline;
        end

        function values_aux = keepSelectedAuxiliary(~, keys_aux, values_aux)
            keep_indices = find(ismember(keys_aux, obj.aux2keep));
            values_aux = values_aux(keep_indices);
        end

        function values_aux = resampleAuxiliary(obj, data, values_aux)
            Fs_resample = round(data.Fs);
            for v = 1:numel(values_aux)
                Fs_i = round(values_aux{v}.Fs);
                values_aux{v}.data = resample(values_aux{v}.data, Fs_resample, Fs_i);
                values_aux{v}.time = data.time;
            end
        end

        function filtered_signal = applyFilter(~, signal, Fs, type, cutoff)
            [b, a] = butter(4, cutoff / (Fs / 2), type);
            filtered_signal = filtfilt(b, a, signal);
        end

        function data = getAuxData(~, keys_aux, values_aux, key)
            index = find(contains(keys_aux, key, 'IgnoreCase', true));
            if ~isempty(index)
                data = values_aux{index}.data;
            else
                data = [];
            end
        end
    end
end
