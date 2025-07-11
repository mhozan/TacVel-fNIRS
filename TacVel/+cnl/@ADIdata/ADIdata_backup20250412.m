classdef ADIdata
    %% ADIdata - Holds ADI Labchart time series data as auxillary data.
    %     data         - Time Series vector(s), e.g. Accelerometer
    %     time         - A columnized vector containing the time
    %     link         - (optional) Linking each column of the data to a string description.
    %     description  - (optional) description of data (e.g. filename)
    %     extra        - (optional) structure containing extra info extracted from ADICHT files:
    %       .downsample_amount
    %       .trigger_minus_rec_start_seconds
    %       .trigger_minus_rec_start_samples
    %       .data_starts
    %       .record_starts
    %       .comments
    %       .units
    %       .data_start_str
    %
    % See also
    % nirs.core.GenericData
    
    properties
        description
        data                        % channel time series in columns
        link                        % object describing geometry
        time                        % vector of time points
        extra                       % struct containing extra info from the ADICHT file
        record_start_delay_s = 0 	% how many seconds after nirx the ADI has started, the first Galileo pulse has triggered ADI recording. Turns out this value is not important.
        stimulus_delay = 0 	        % how many seconds after nirx recoding started, the first Galileo pulse has arrived.
        discardthedelay = false    	% true if the time and data vectors are to be synchronized with the nirx data
        resampled = false           % true if the time and data vectors are resampled post-recording.
        keepit = true               % true if the entire aux signal is to be kept, false if to be ignored/discarded.
      end
    
    properties( Dependent = true )
        Fs                          % sampling frequency in Hz
    end

    methods
        function obj = ADIdata( data, time,link, description , extra , record_start_delay_s , stimulus_delay, discardthedelay )
            %% Data - Creates a Data object.
            % 
            % Args:
            %     data         - Time Series vector(s), e.g. Accelerometer
            %     time         - A columnized vector containing the time
            %     link         - (optional) Linking each column of the data to a string description.
            %     description  - (optional) description of data (e.g. filename)
            %     extra        - (optional) structure containing extra info extracted from ADICHT files:
            %       .downsample_amount
            %       .trigger_minus_rec_start_seconds
            %       .trigger_minus_rec_start_samples
            %       .data_starts
            %       .record_starts
            %       .comments         
            %       .units            
            %       .data_start_str

            

            if nargin > 0, obj.data         = data;         end
            if nargin > 1, obj.time         = time;         end
            if nargin > 2, obj.link        = link;        end
            if nargin > 3, obj.description  = description;  end
            if nargin > 4, obj.extra  = extra;  end
            if nargin > 5, obj.record_start_delay_s  = record_start_delay_s;  end
            if nargin > 6, obj.stimulus_delay  = stimulus_delay;  end
            if nargin > 7, obj.discardthedelay  = discardthedelay;  end
            
            %% Process some auxiliary data
            auxtobeprocessed = {'mic';'eyelid';'hand';'ADI';'accl';'Resp';'Galileo'};
            auxtobeprocessed = {'eyelid';'accl'};
            auxtobeprocessed = {'eyelid'};
            strmtch=regexpi(obj.description,auxtobeprocessed);
            idx = find(~cellfun(@isempty,strmtch)); 
            if contains(obj.description,'eyelid')
                        eyeblink = emg2blinkrate(obj.data);
                        jawmovement = (highpass(obj.data,20,obj.Fs,'ImpulseResponse','fir'));
                        data = [eyeblink,jawmovement];
                        obj.data=data;
            elseif contains(obj.description,'accl')
                % obj.data = lowpass(obj.data,200,obj.Fs,'ImpulseResponse','fir');
                % obj.data = sqrt(highpass(obj.data,0.1,obj.Fs,'ImpulseResponse','fir').^2);

                [b, a] = butter(2, [.5 20] / (obj.Fs / 2), 'bandpass');
                obj.data = filtfilt(b, a, obj.data);
                % figure, 
                % plot(obj.data)
                % hold on
                % obj.data = sqrt(obj.data.^2);
                % plot(obj.data,'--r')

            elseif contains(obj.description,'mic')
                        mic = (highpass(obj.data,100,obj.Fs,'ImpulseResponse','fir'));
                        mic_rms = sqrt(mic.^2);
                        obj.data = mic_rms;%movmean(mic_rms,1*obj.Fs);
            elseif contains(obj.description,'hand')
                        handmovement = (highpass(obj.data,20,obj.Fs,'ImpulseResponse','fir'));
                        handmovement = sqrt(handmovement.^2);
                        handmovement = movmean(handmovement,obj.Fs);
                        obj.data = handmovement;
            elseif contains(obj.description,'ADI')
                        obj.data = heartbeatestimate(obj.data);

            elseif contains(obj.description,'Resp')
                        [peaks,~,~] = peaks2rate(obj.data, ...
                            'Fs',obj.Fs,...
                            'f_lowpass',5, ...
                            'MinPeakProminence',.2, ...
                            'MinPeakDistance',3*obj.Fs, ... %max 20 breaths per minute
                            'MinPeakWidth',1.5*obj.Fs, ...
                            'MinPeakHeight',-inf);
                        % figure(11), clf; plot(obj.time,obj.data), hold on, scatter(y2(y3),ones(size(y3))),plot(obj.time,60*peaks)
                        obj.data = movmean(60*peaks,60*obj.Fs);
            % elseif contains(obj.description,'Galileo')
            else
                %do nothing

            end


            if ~isempty(idx) && false
                % obj.data=zscore(obj.data);
                % figure(11), clf; plot(obj.time,obj.data), hold on, 
                switch idx
                    case 1 %mic
                        mic = (highpass(obj.data,100,obj.Fs,'ImpulseResponse','fir'));
                        mic_rms = sqrt(mic.^2);
                        obj.data = mic_rms;%movmean(mic_rms,1*obj.Fs);

                    case 2 %EMG_eyelid
                        % disp(obj.extra.ADICHTfilepath)
                        eyeblink = emg2blinkrate(obj.data);
                        % lowpass(EMG_eyelid,f_lowpass,Fs,'ImpulseResponse','fir')
                        jawmovement = (highpass(obj.data,20,obj.Fs,'ImpulseResponse','fir'));
                        % bhi = fir1(2^10,[0.1 ]*2/obj.Fs);
                        % jawmovement = jawmovement - mean(jawmovement);
                        % jawmovement = sqrt(jawmovement.^2);
                        % jawmovement = detrend_lin_pw(jawmovement,round(100*obj.Fs));
                                            % d = detrend_lin_pw(d,round(obj.detrendpiecewiseduration*data(i).Fs));


                    % d = filtfilt(bhi,1,d);

                        % jawmovement = sqrt(jawmovement.^2);
                        % jawmovement = movmean(jawmovement,obj.Fs);
                        data = [eyeblink,jawmovement];
                        % [eyeblink,y2,y3] = emg2blinkrate(obj.data);
                        % figure(11), clf; plot(obj.time,obj.data), hold on, scatter(y2(y3),ones(size(y3))),plot(obj.time,eyeblink)
                        obj.data=data;

                    case 3 %EMG_Hand
                        handmovement = (highpass(obj.data,20,obj.Fs,'ImpulseResponse','fir'));
                        handmovement = sqrt(handmovement.^2);
                        handmovement = movmean(handmovement,obj.Fs);
                        obj.data = handmovement;
                       


                    case 4 %ADI Pulse (Heart rate)
                        obj.data = heartbeatestimate(obj.data);
                        obj.data = obj.data/60; %convert to beats per seconds

                    case 5 %Acceleromter
                        obj.data = lowpass(obj.data,50,obj.Fs,'ImpulseResponse','fir');
                        obj.data = sqrt(highpass(obj.data,1,obj.Fs,'ImpulseResponse','fir').^2);

                        % obj.data = lowpass_accelerometer(obj.data,obj.Fs);

                    case 6 %Respiration
                        [peaks,~,~] = peaks2rate(obj.data, ...
                            'Fs',obj.Fs,...
                            'f_lowpass',5, ...
                            'MinPeakProminence',.2, ...
                            'MinPeakDistance',3*obj.Fs, ... %max 20 breaths per minute
                            'MinPeakWidth',1.5*obj.Fs, ...
                            'MinPeakHeight',-inf);
                        % figure(11), clf; plot(obj.time,obj.data), hold on, scatter(y2(y3),ones(size(y3))),plot(obj.time,60*peaks)
                        obj.data = movmean(60*peaks,60*obj.Fs);
                        % obj.data = (highpass(obj.data,.01,obj.Fs));
                        % plot(obj.time,ata)

                    case 7 %galielo
                        % obj.data = 1*obj.data;


                    otherwise
                        keyboard

                end
            end
%%
            %
            %            
            if obj.discardthedelay
                obj.time = obj.time+stimulus_delay+extra.trigger_minus_rec_start_seconds;

                % obj.time = obj.time+record_start_delay_s+extra.trigger_minus_rec_start_seconds;
%                 [~,discardindex] = min(abs(obj.time - obj.record_start_delay_s));
%                 obj.time(1:discardindex-1) = [];
%                 obj.time = obj.time-obj.time(1);
%                 obj.data(1:discardindex-1,:) = [];
                
            end

        end
        
        
        function obj = set.time( obj, time )
           assert( isvector(time) )
           obj.time = time(:);
        end
        
     	function vec = getStimVector( obj, time )
            vec = interp1( obj.time, obj.data, time );
        end
        
        function out = get.Fs( obj )
            if length(obj.time) > 1
                out = 1 / mean(diff( obj.time ));
            else
                out = NaN;
            end
        end
                
        function out = sorted( obj, colsToSortBy )
            %% returns sorted channels of data by column in probe.link
            out = obj;
            if nargin < 2
                colsToSortBy = {'name', 'type'};
            end
            if(length(obj)>1)
                for idx=1:length(obj)
                    out(idx)=sorted(obj(idx),colsToSortBy );
                end
                return
            end
            [out.link, idx] = sortrows(out.link, colsToSortBy);
            out.data = out.data(:,idx);
        end
        
        function varargout=draw( obj, lstChannels,adderr )
            %% draw - Plots the probe geometry.
            % 
            % Args:
            %     lstChannels - list of channels to show
            
            % get data
            if nargin == 1
                lstChannels = 1:size(obj(1).data,2);
                
            end
            if(nargin<3)
                adderr=false;
            end
            
            if(isempty(lstChannels))
                % draw a plot, but then turn it off to only show the stim
                % marks
                lstChannels=1;
                showplot=false;
            else
                showplot=true;
            end
            if(isempty(obj(1).data))
                lstChannels=[];
            end
            
            if(length(obj)>1)
                figure;
                a=round(sqrt(length(obj)));
                b=ceil(sqrt(length(obj)));
                for i=1:length(obj)
                    subplot(a,b,i);
                    obj(i).draw(lstChannels,adderr);
                    legend off;
                end
                return
            end
            
            t = obj.time;
            d = obj.data(:,lstChannels);
            
            % plots
            gca; hold on;
            
            % data min/max/size
            dmax = max( real(d(:)) );
            dmin = min( real(d(:)) );
            dsize = (dmax-dmin);
            
          
                % min/max of axes
             	pmin = dmin - 0.1*dsize;
                pmax = dmax + 0.1*dsize;
            
            % plot data
            if(~isreal(d) & adderr)
                h=errorbar( t, real(d),imag(d) );
            else
                h=plot( t, real(d) );
            end
            xlabel( 'seconds' );

            
            % axes limits
            xlim( [min(t) max(t)] );
            if pmin == pmax
                ylim(pmin + [-1 1]);
            else
                ylim( [pmin pmax] );
            end
            
            if(~showplot)
                set(h,'visible','off');
            end
            
            if(nargout>0)
                varargout{1}=h;
            end
            
        end
        
    end
end

