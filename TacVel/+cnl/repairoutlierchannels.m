classdef repairoutlierchannels < nirs.modules.AbstractModule
    %% repairoutlierchannels applies a 0.02-0.25 Hz BPF only to outlier channels.
    %
    %   Identifies outlier channels based on high variability (moving std dev + MAD)
    %   and applies a specific bandpass filter (defined by hipasscutoff_Hz and
    %   lowpasscutoff_Hz) ONLY to these channels using matrix operations.
    %   Throws an error if data length is too short for the specified FIR filter order.
    %
    %   Properties:
    %       deviation_threshold (double): Factor for MAD threshold to detect outliers. Default: 3.
    %       outlier_detection_window_sec (double): Window duration (sec) for moving std dev. Default: 40.
    %       hipasscutoff_Hz (double): Lower cutoff for the BPF applied *only* to outliers. Default: 0.01.
    %       lowpasscutoff_Hz (double): Upper cutoff for the BPF applied *only* to outliers. Default: 0.025.
    %       fir1ord (integer): Order for the FIR bandpass filter. Default: 2048. Must satisfy N > 3*order.
    %       keep_mean (logical): If true, add original channel mean back after filtering. Default: false.

    properties
        deviation_threshold         = 10;     % Threshold for MAD of moving std dev
        outlier_detection_window_sec= 20;    % Window duration for movstd in seconds
        hipasscutoff_Hz             = 0.02;  % BPF lower cutoff [Hz] for outliers
        lowpasscutoff_Hz            = 0.25; % BPF upper cutoff [Hz] for outliers
        fir1ord                     = 4095;  % FIR filter order
        keep_mean                   = true; % Add original DC offset back after filtering?
        dovisualize                 = true;
    end

    methods
        function obj = repairoutlierchannels( prevJob )
            obj.name = 'Repair Outlier Channels';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end

        function data = runThis( obj, data )
            for i = 1:numel(data)
                d_raw = data(i).data;
                Fs = data(i).Fs;
                n_samples = size(d_raw, 1);
                n_channels = size(d_raw, 2);

                % fprintf('Processing data entry %d/%d...\n', i, numel(data));

                % --- Outlier Detection ---
                windowSizeSamples = round(obj.outlier_detection_window_sec * Fs);
                windowSizeSamples = max(2, min(windowSizeSamples, n_samples - 1));
                movingStd = movstd(d_raw, windowSizeSamples);
                movingStd = movingStd;
                stdevValues = std(movingStd);%, 1, 1, 'omitnan');
                % rmsValues = rms(movingStd);
                medianSTDEV = median(stdevValues);
                % medianRMS = median(rmsValues);
                isOutlierChannel = stdevValues > (obj.deviation_threshold * medianSTDEV);
                outlier_indices = find(isOutlierChannel);

                figure(10001),clf, subplot(4,1,1), plot(d_raw),
                subplot(4,1,2), plot(movingStd(:,~isOutlierChannel)), 
                subplot(4,1,3), plot(movingStd(:,outlier_indices)), 
                subplot(4,1,4), scatter(1:length(stdevValues),stdevValues), hold on, plot(1:122,obj.deviation_threshold*medianSTDEV*ones(1,122))
                % disp(nnz(isOutlierChannel))
                drawnow


                % --- Selective Filtering (Vectorized) ---
                if ~isempty(outlier_indices)
                    fprintf('  Found %d outlier channels: %s\n', length(outlier_indices), num2str(outlier_indices'));

                    % --- Check Data Length vs FIR Order ---
                    % Removed the while loop for dynamic order reduction.
                    % Now, it errors if data is too short for the specified order.
                    min_data_len = 3 * obj.fir1ord + 1; % filtfilt requirement
                    if n_samples < min_data_len
                        error('Data length (%d samples) is too short for the specified FIR filter order (%d). Minimum required is %d (3*order+1). Choose a smaller fir1ord or use longer data segments.', n_samples, obj.fir1ord, min_data_len);
                    end

                    % --- Design BPF ---
                    nyquist = Fs / 2;
                    low_cutoff = obj.hipasscutoff_Hz / nyquist;
                    high_cutoff = obj.lowpasscutoff_Hz / nyquist;
                    if low_cutoff <= 0 || high_cutoff >= 1 || low_cutoff >= high_cutoff
                        error('Invalid BPF cutoff frequencies [%f, %f] Hz for Fs = %f Hz.', obj.hipasscutoff_Hz, obj.lowpasscutoff_Hz, Fs);
                    end
                    bpf_coeffs = fir1(obj.fir1ord, [low_cutoff, high_cutoff], 'bandpass');

                    % --- Apply Filter (Vectorized) ---
                    outlier_data = d_raw(:, outlier_indices);
                    original_means = [];

                    if obj.keep_mean
                        original_means = mean(outlier_data, 1, 'omitnan');
                    end


                    filtered_outliers = filtfilt(bpf_coeffs, 1, outlier_data);

                    if obj.keep_mean && ~isempty(original_means)
                        filtered_outliers = filtered_outliers + original_means;
                    end

                    d_processed = d_raw;
                    d_processed(:, outlier_indices) = filtered_outliers;
                    data(i).data = d_processed; % Update data object *only if successful*
                    fprintf('  Successfully applied BPF to outlier channels.\n');

                    % zzz= find(data(i).data<0);
                    if any(data(i).data(:)<0)
                        % keyboard
                        i
                        absoluteminimumvalue = min(data(i).data(:))
                        data(i).data = data(i).data - absoluteminimumvalue + eps;
                    end
                    % data

                    % ... inside the 'if ~isempty(outlier_indices)' block ...

                    % --- Apply Filter (Vectorized) ---
                    % (try-catch block for filtering happens here)
                    % At this point, 'd_processed' holds the potentially modified data.
                    % If filtering failed, 'd_processed' might still just be 'd_raw'.

                    % --- Visualization ---
                    if obj.dovisualize && exist('d_processed', 'var') % Check if d_processed exists
                        figure(900 + i); % Unique figure per data entry
                        clf;
                        t = data(i).time;
                        num_channels = size(d_raw, 2);
                        channel_indices = 1:num_channels;
                        non_outlier_indices = setdiff(channel_indices, outlier_indices); % More robust way

                        % Get threshold (recalculate if meanMad was NaN/zero earlier)
                        threshold_val = NaN;
                        if isfinite(medianSTDEV) && medianSTDEV > 0
                            threshold_val = obj.deviation_threshold * medianSTDEV;
                        end

                        % 1. MAD values vs Threshold
                        ax1 = subplot(4, 1, 1);
                        scatter(ax1, channel_indices, stdevValues, 'b.'); hold on;
                        scatter(ax1, outlier_indices, stdevValues(outlier_indices), 'r.');
                        if isfinite(threshold_val)
                            plot(ax1, xlim, [threshold_val, threshold_val], 'r--', 'LineWidth', 1);
                        end
                        hold off; grid on;
                        title(ax1, 'MAD of Moving Std Dev vs. Channel');
                        ylabel(ax1, 'MAD Value');
                        legend(ax1, 'Non-Outlier', 'Outlier', 'Threshold', 'Location', 'best');
                        set(ax1, 'XTickLabel', []); % Remove x-ticks for this plot

                        % 2. Non-Outlier Channels (Raw)
                        ax2 = subplot(4, 1, 2);
                        if ~isempty(non_outlier_indices)
                            plot(ax2, t, d_raw(:, non_outlier_indices));
                            ylabel(ax2, 'Amplitude'); title(ax2, 'Non-Outlier Channels (Raw)'); grid on;
                        else
                            title(ax2, 'No Non-Outlier Channels'); axis off; % Handle case if all are outliers
                        end
                        set(ax2, 'XTickLabel', []); % Remove x-ticks for this plot


                        % 3. Outlier Channels (Raw)
                        ax3 = subplot(4, 1, 3);
                        plot(ax3, t, d_raw(:, outlier_indices));
                        ylabel(ax3, 'Amplitude'); title(ax3, [num2str(length(outlier_indices)),' Outlier Channels (Raw)']); grid on;
                        set(ax3, 'XTickLabel', []); % Remove x-ticks for this plot


                        % 4. Outlier Channels (Processed)
                        ax4 = subplot(4, 1, 4);
                        plot(ax4, t, d_processed(:, outlier_indices)); % Plot from d_processed
                        ylabel(ax4, 'Amplitude'); title(ax4, 'Outlier Channels (Filtered/Processed)'); grid on;
                        xlabel(ax4, 'Time (s)');

                        % Link time axes
                        linkaxes([ax2, ax3, ax4], 'x');
                        if ~isempty(t)
                            xlim(ax2,[t(1) t(end)]); % Set x limits based on time
                        end

                        % Super Title for the figure
                        sgtitle(sprintf('Data Entry %d: Outlier Repair Visualization (Threshold Factor: %.2f)', i, obj.deviation_threshold), 'Interpreter', 'none');

                        drawnow; % Update figure window
                    end

                    % ... rest of the code for loop i ...
                else
                    fprintf('  No outlier channels detected.\n');
                    % No changes needed
                end
            end
            fprintf('Finished processing all data entries.\n');
        end
    end
end