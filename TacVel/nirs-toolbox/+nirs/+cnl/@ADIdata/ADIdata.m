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
            if contains(obj.description,'eyelid')
            fh= 0.5; fl =10; filter_order= 1023; 
                        eyeblink_rate = emg2blinkrate(obj.data);
                        eyeblink = applyBandPassFilter(obj.data,...
                            'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);
            fh= 10; fl =100; filter_order= 1023; 
                        facialEMG = applyBandPassFilter(obj.data,...
                            'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);
                        % facialEMG = (highpass(obj.data,20,obj.Fs,'ImpulseResponse','fir'));
                        eyeface = [eyeblink,facialEMG,eyeblink_rate];
                        % figur(obj.description), clf, plot(time,eyeface), hold on, plot(time,obj.data)
                        obj.data=eyeface;
            elseif contains(obj.description,'accl')
                % obj.data = lowpass(obj.data,200,obj.Fs,'ImpulseResponse','fir');
                % obj.data = sqrt(highpass(obj.data,0.1,obj.Fs,'ImpulseResponse','fir').^2);
                fh= 0.5; fl =10; filter_order= 4095; 
                accl_f = applyBandPassFilter(obj.data, 'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);
                % figur(obj.description), clf, plot(time,accl_f), hold on, plot(time,obj.data-1.5)

                    %% x
                % [b, a] = butter(2, [.5 20] / (obj.Fs / 2), 'bandpass');
                % obj.data = filtfilt(b, a, obj.data);
                obj.data = accl_f;
                % figure, 
                % plot(obj.data)
                % hold on
                % obj.data = sqrt(obj.data.^2);
                % plot(obj.data,'--r')

            elseif contains(obj.description,'mic')
                        fh= 2; fl =10; filter_order= 1023; 
                        mic = applyBandPassFilter(obj.data, 'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);
                        % figur(obj.description), clf, plot(time,mic-100), hold on, plot(time,obj.data)
                        obj.data = mic;%movmean(mic_rms,1*obj.Fs);
            elseif contains(obj.description,'hand')
                        fh= 0.5; fl =10; filter_order= 1023;
                        handmovement = applyBandPassFilter(obj.data,...
                            'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);                       
                % handmovement = (highpass(obj.data,20,obj.Fs,'ImpulseResponse','fir'));
                % handmovement = sqrt(handmovement.^2);
                % handmovement = movmean(handmovement,obj.Fs);
                        % figur(obj.description), clf, plot(time,obj.data), hold on, plot(time,handmovement)
                        obj.data = handmovement;
            elseif contains(obj.description,'ADI')
                        % obj.data = heartbeatestimate(obj.data);
                        fh= 0.5; fl =4; filter_order= 1023; 
                        obj.data = applyBandPassFilter(obj.data,...
                            'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);
                        % figur(obj.description), clf, plot(time,data), hold on, plot(time,ata)
            elseif contains(obj.description,'Resp')
                        % [peaks,~,~] = peaks2rate(obj.data,'Fs',obj.Fs,'f_lowpass',5,'MinPeakProminence',.2,'MinPeakDistance',3*obj.Fs,'MinPeakWidth',1.5*obj.Fs,'MinPeakHeight',-inf);
                        fh= .1; fl =2; filter_order= 1023;
                        resp = applyBandPassFilter(obj.data,...
                            'lowpass', fl, 'highpass', fh, 'Fs', obj.Fs,'filter_order',filter_order);
                        % figur(obj.description),clf; plot(obj.time,obj.data), hold on,plot(obj.time,resp), title(obj.description)
                        % plot(obj.time,resp_detrend+2)%,plot(obj.time,movmean(60*peaks,60*obj.Fs)), 
                        % obj.data = movmean(60*peaks,60*obj.Fs);
                        obj.data = resp;
            % elseif contains(obj.description,'Galileo')
            else
                %do nothing

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

