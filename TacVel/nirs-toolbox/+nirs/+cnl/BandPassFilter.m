classdef BandPassFilter < nirs.modules.AbstractModule
%% Simple Band Pass Filter (FIR + filtfilt)
% Options:
%   - lowpass: upper frequency (Hz)
%   - highpass: lower frequency (Hz)
%   - do_downsample: optionally resample based on lowpass
%   - keepdc: retain mean of signal if true

    properties
        lowpass;
        highpass;
        do_downsample;
        keepdc;
        filter_order = 510; % FIR filter order
        makeitpositive; %make the output ready for log fn in nirs.modules.OpticalDensity
    end

    methods
        function obj = BandPassFilter(prevJob)
            obj.name = 'FIR Band Pass Filter nirs_cnl_bandpassfilter';
            obj.lowpass = [];
            obj.highpass = [];
            obj.do_downsample = false;
            obj.keepdc = true;
            obj.makeitpositive =false;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end

        function data = runThis(obj, data)
            for i = 1:length(data)

                Fs = data(i).Fs;

                % Check cutoff presence
                if isempty(obj.lowpass) && isempty(obj.highpass)
                    return;
                end

                % Normalize frequency cutoff
                Wn = [];
                if ~isempty(obj.highpass)
                    Wn(1) = obj.highpass / (Fs/2);
                end
                if ~isempty(obj.lowpass)
                    Wn(end+1) = obj.lowpass / (Fs/2);
                end

                % Design FIR filter
                if length(Wn) == 1
                    if isempty(obj.highpass)
                        b = fir1(obj.filter_order, Wn, 'low');
                    else
                        b = fir1(obj.filter_order, Wn, 'high');
                    end
                else
                    b = fir1(obj.filter_order, Wn, 'bandpass');
                end
                % figur(mfilename),clf, freqz(b,1,obj.filter_order)

                % Apply filtering
                d = data(i).data;

                % Remove DC
                if obj.keepdc
                    dc = mean(d,1);
                else
                    dc = zeros(1, size(d,2));
                end
                d = d - dc;

                % Zero-phase filtering
                d = filtfilt(b, 1, d);

                % Restore DC if needed
                if obj.keepdc
                    d = d + dc;
                end
                if obj.makeitpositive
                    d = d+1+min(d(:));
                end

                % Optional downsampling
                if obj.do_downsample && ~isempty(obj.lowpass) && obj.lowpass*2 < Fs
                    t = data(i).time;
                    N = floor((t(end)-t(1)) * obj.lowpass * 2);
                    new_t = t(1) + (0:N-1)' / (obj.lowpass * 2);
                    d = interp1(t, d, new_t, 'linear', 'extrap');
                    data(i).time = new_t;
                end

                % Replace data
                data(i).data = d;
            end
        end
    end
end
