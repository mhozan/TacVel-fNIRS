classdef pruneaux < nirs.modules.AbstractModule
    %% PRUNEAUX prunes the aux channels
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
        % what names to be kept as aux regressors, rest will be thrown out.
        aux2keep = {'ADIpulsetransducer',...heart rate
            'mic',...
            'Eyeblink',...
            'iblink_rate',...
            'EMG_face',...
            'hand',... %EMG2_hand or EMG2_righthand
            ...'XYZ_acceleration',... x+y+z or x*y*z or another combination of sort
            'Xaccel',...
            'Yaccel',...
            'Zaccel',...
            'Respiration_rate'};%,...
            % 'CombinedAuxMotion',...
            % 'NIRmidline'};
        % resampleaux=false;
        resampleaux=false;
        trimtimevectors = true;
        clip_threshold = 99;
        visualize=false;
    end
    
    methods
        function obj = pruneaux( prevJob )
            obj.name = 'nirs_cnl_pruneaux';

            if nargin > 0
                obj.prevJob = prevJob;
            end
        end

        function data = runThis( obj, data )
            indexaux2keep = zeros(size(obj.aux2keep));
            for i = 1:numel(data)
                if isempty(data(i).auxillary)
                    return
                end

                keys_aux = data(i).auxillary.keys;
                values_aux = data(i).auxillary.values;

                %% accelerometer
                acclindices = find(contains(keys_aux,'accl','IgnoreCase',true));%this should have 3 or 0 matches
                if ~isempty(acclindices)
                    combinedaccelerometer = values_aux(acclindices); 
                   
                    %% Bandpass filter specs 
                    fir1order=511;
                    fl=10;
                    fh=2;
                    Fs = combinedaccelerometer{1}.Fs; %should be the same for all 3 accel channels

                    %% x
                    x = combinedaccelerometer{1}.data(:);                 % x = sqrt(x.^2);
                    % x_f = applyBandPassFilter(x, 'lowpass', fl, 'highpass', fh, 'Fs', Fs,'filter_order',fir1order);
                    values_aux(end+1) = values_aux(acclindices(1)); %duplicate accl_X
                    keys_aux(end+1) = {'Xaccel'}; %change its name
                    values_aux{end}.description = 'Xaccel';
                    values_aux{end}.data = x;

                    %% y 
                    y=combinedaccelerometer{2}.data(:);
                    % y_f = applyBandPassFilter(y, 'lowpass', fl, 'highpass', fh, 'Fs', Fs,'filter_order',fir1order);
                    values_aux(end+1) = values_aux(acclindices(2)); %duplicate accl_Y
                    keys_aux(end+1) = {'Yaccel'}; %change its name
                    values_aux{end}.description = 'Yaccel';
                    values_aux{end}.data = y;

                    %% z
                    z=combinedaccelerometer{3}.data(:);
                    % z_f = applyBandPassFilter(z, 'lowpass', fl, 'highpass', fh, 'Fs', Fs,'filter_order',fir1order);
                    values_aux(end+1) = values_aux(acclindices(3)); %duplicate accl_Z
                    keys_aux(end+1) = {'Zaccel'}; %change its name
                    values_aux{end}.description = 'Zaccel';
                    values_aux{end}.data = z;

                    %% x+y+z
                    % xyz_f = 1000*applyBandPassFilter(x.*y.*z, 'lowpass', fl, 'highpass', fh, 'Fs', Fs,'filter_order',fir1order);
                    xyz_f = applyBandPassFilter(x+y+z, 'lowpass', fl, 'highpass', fh, 'Fs', Fs,'filter_order',fir1order);
                    % xyz_f = applyBandPassFilter(sqrt(x.^2+y.^2+z.^2), 'lowpass', fl, 'highpass', fh, 'Fs', Fs,'filter_order',fir1order);
                    % xyz_f = sqrt(x_f.^2+y_f.^2+z_f.^2);
                    values_aux(end+1) = values_aux(acclindices(1)); %duplicate accl_X
                    keys_aux(end+1) = {'XYZ_acceleration'}; %change its name
                    values_aux{end}.description = 'XYZ_acceleration';
                    values_aux{end}.data = xyz_f;
                    
                    % figur(['xyz_accl_',num2str(i,'%0.2d')]), clf
                    % tvec = values_aux{end}.time;
                    % l1 = plot(tvec, x     +.01 ); l1.DisplayName = 'x_raw_acc'; hold on
                    % l2 = plot(tvec, y     +.11); l2.DisplayName = 'y_raw_acc';
                    % l3 = plot(tvec, z     +.21); l3.DisplayName = 'z_raw_acc';
                    % l4 = plot(tvec, x+y+z +.31); l4.DisplayName = 'x+y+z_raw_acc'; 
                    % l5 = plot(tvec, x_f    ); l5.DisplayName = 'x_f_acc'; hold on
                    % l6 = plot(tvec, y_f   +.1); l6.DisplayName = 'y_f_acc';
                    % l7 = plot(tvec, z_f   +.2); l7.DisplayName = 'z_f_acc';
                    % l8 = plot(tvec, xyz_f +.3); l8.DisplayName = 'x+y+z_f_acc'; 
                    % legend('Interpreter','none')
                end                

                %% Respiration
                % respindices = find(contains(keys_aux,{'resp','cage','abdomen'},'IgnoreCase',true));
                respindices = find(contains(keys_aux,{'_rclo','ribcagelo'},'IgnoreCase',true));
                if ~isempty(respindices)
                    RESP_aux = values_aux{respindices}.data;
                    % RESP_aux_combined = [RESP_aux{1}.data,...
                    %     RESP_aux{2}.data,RESP_aux{3}.data,...
                    %     RESP_aux{4}.data];
                    % normalizeeachrespcolumn = true;
                    % if normalizeeachrespcolumn
                    %     respnormalizationfactor  = prctile(RESP_aux_combined,90);
                    %     if all(respnormalizationfactor) % all resp channels are mostly non zero
                    %     RESP_aux_combined = RESP_aux_combined./respnormalizationfactor;
                    %     else
                    %     RESP_aux_combined = 0*RESP_aux_combined;
                    %     disp(['i = ',num2str(i,'%0.2d'),'Respiration rate data is not good. set to zero.'])
                    %     end
                    % end
                    % 
                    % RESP_median = median(RESP_aux_combined,2)./60;
                    % 

                    % figur(['respiration_',num2str(i,'%0.2d')]), clf
                    % plot(RESP_aux_combined), hold on
                    % plot(60*RESP_median)
                    % legend
                    % RESP_median = lowpass(RESP_median,20,values_aux{respindices(1)}.Fs,'ImpulseResponse','fir');
                    % RESP_median = highpass(RESP_median,.1,values_aux{respindices(1)}.Fs,'ImpulseResponse','fir');

                    values_aux(end+1) = values_aux(respindices(1)); %duplicate RCLo_Resp
                    keys_aux(end+1) = {'Respiration_rate'}; %change its name
                    values_aux{end}.description = 'Respiration_rate';
                    values_aux{end}.data = RESP_aux; %replace its data

                    % RESP_median(RESP_median>prctile(RESP_median,obj.clip_threshold)) = prctile(RESP_median,obj.clip_threshold);
                    % auxmultiplied =  auxmultiplied.*RESP_median/(max(RESP_median));
                end

                %% Separate eyeblink and facial muscle tone and eyeblink rate
                eyeEMGindex = find(contains(keys_aux,'eye','IgnoreCase',true));
                if ~isempty(eyeEMGindex)
                    eyeEMGdata = values_aux{eyeEMGindex}.data;
                    eyeblink = eyeEMGdata(:,1);
                    facialEMG = (eyeEMGdata(:,2).^2);
                    eyeblink_rpm = eyeEMGdata(:,3);
                    % facialEMG = eyeEMGdata(:,2);

                    values_aux(end+1) = values_aux(eyeEMGindex); %duplicate 'CH05_07_EMG1_eyelid'
                    keys_aux(end+1) = {'Eyeblink'}; %change its name
                    values_aux{end}.description = 'Eyeblink';
                    values_aux{end}.data = eyeblink; %replace its data

                    values_aux(end+1) = values_aux(eyeEMGindex); %duplicate 'CH05_07_EMG1_eyelid'
                    keys_aux(end+1) = {'EMG_face'}; %change its name
                    values_aux{end}.description = 'EMG_face';
                    % facialEMG2= facialEMG;
                    % facialEMG2((round(.4*end):round(.6*end)))=0;
                    % facialEMG2((1:round(.05*end)))=0;
                    % facialEMG2(round(.95*end):end)=0;
                    % threshemg = prctile(facialEMG2,obj.clip_threshold);
                    % facialEMG(facialEMG2>threshemg) = threshemg;
                    % facialEMG = facialEMG./threshemg;
                    % facialEMG = highpass(facialEMG,5,values_aux{eyeEMGindex}.Fs,'ImpulseResponse','fir');
                    % facialEMG = facialEMG.^2;
                    % 
                    values_aux{end}.data = facialEMG; %replace its data
                    % auxmultiplied =  auxmultiplied.*facialEMG/(max(facialEMG));

                    values_aux(end+1) = values_aux(eyeEMGindex); %duplicate 'CH05_07_EMG1_eyelid'
                    keys_aux(end+1) = {'iblink_rate'}; %change its name
                    values_aux{end}.description = 'iblink_rate';
                    values_aux{end}.data = eyeblink_rpm; %can be used for detecting sleep

                else
                    fprintf('no eye_emg auxiliary data found.')
                end
                %% rename EMG_hand
                emghandindex = find(contains(keys_aux,'hand','IgnoreCase',true));
                if ~isempty(emghandindex)
                    values_aux{emghandindex}.description = 'Hand_EMG2';
                    keys_aux{emghandindex} = 'Hand_EMG2';
                end
                %% rename ADI
                adi_index = find(contains(keys_aux,'adi','IgnoreCase',true));
                if ~isempty(adi_index)
                    values_aux{adi_index}.description = 'Heartrate_ADI_pulse';
                    keys_aux{adi_index} = 'Heartrate_ADI_pulse';
                end
                %% CombinedAuxMotion: a multiplication of a few auxiliary channels to best capture motion artifacts in a single channel
                % % clip_threshold = 99;
                %  emghandindex = find(contains(keys_aux,'hand','IgnoreCase',true));
                %  emghand = values_aux{emghandindex}.data;
                %  emghand(emghand>prctile(emghand,obj.clip_threshold)) = prctile(emghand,obj.clip_threshold);
                %  % auxmultiplied = auxmultiplied.*emghand/(max(emghand));
                % 
                %  micindex = find(contains(keys_aux,'mic','IgnoreCase',true));
                %  micdata = values_aux{micindex}.data;
                %  micdata = highpass(micdata,10,values_aux{micindex}.Fs,'ImpulseResponse','fir');
                %  micdata(micdata>prctile(micdata,obj.clip_threshold)) = prctile(micdata,obj.clip_threshold);
                %  % auxmultiplied = auxmultiplied.*micdata/(max(micdata));
                % 
                %  % thresh = prctile(auxmultiplied,obj.clip_threshold);
                %  % auxmultiplied(auxmultiplied>thresh)=thresh;
                %  % auxmultiplied=auxmultiplied./thresh;
                %  fc_hipass= 1.5;
                %  % auxmultiplied = highpass(auxmultiplied,fc_hipass,values_aux{1}.Fs,'ImpulseResponse','fir');
                % 
                % 
                %  values_aux(end+1) = values_aux(1); %duplicate an aux
                %  keys_aux(end+1) = {'CombinedAuxMotion'}; %change its name
                %  values_aux{end}.description = 'CombinedAuxMotion';
                %  % values_aux{end}.data = auxmultiplied; %replace its data'
                %  % combinedmotion =  sqrt(0.5*(facialEMG+xyzaccl).^2); %replace its data;
                %  % combinedmotion =  facialEMG+xyz_acc_env+micdata+emghand; %replace its data;
                %  combinedmotion =  facialEMG+XplusYplusZ_filtered+micdata+emghand; %replace its data;
                %  % combinedmotion = highpass(combinedmotion,.5,values_aux{1}.Fs,'ImpulseResponse','fir');
                %  combinedmotion = combinedmotion.^2;
                %  combinedmotion = highpass(combinedmotion,10,values_aux{1}.Fs,'ImpulseResponse','fir');
                %  combinedmotion = abs(combinedmotion);
                %  combinedmotion = combinedmotion - min(combinedmotion);
                %  combinedmotion((combinedmotion>prctile(combinedmotion,90)))  = prctile(combinedmotion,90);
                % 
                %  values_aux{end}.data = combinedmotion; %replace its data

                 %% Consider centered optode links (S8-D9, S9-D9, S9-D10) as regressors
                 % 
                 % lst = ismember(data(i).probe.link.source,[8,9])  & ismember(data(i).probe.link.detector,[9,10]);
                 % if isempty(lst) || ~nnz(lst)
                 %     keyboard
                 % end
                 % if ~isempty(lst) && nnz(lst)
                 %     midline = data(i).data(:,lst);
                 %     midline1= midline(:,1);
                 %     midline2= midline(:,2);
                 %     midline3= midline(:,3);
                 %     midline4= midline(:,4);
                 %     midline5= midline(:,5);
                 %     midline6= midline(:,6);
                 %     % midline = mean(midline.').'; %unstable model?
                 %     midline = median(midline.').';
                 %     midline = midline-min(midline)+eps; %make sure deta is greater than zero
                 % 
                 %     %optical density
                 %     m = mean( midline, 1 ,"omitnan");
                 %     midline =log(midline./(ones(size(midline,1),1)*m));
                 %     midline(midline<0)=0;
                 %     if any(isnan(midline)) || any(isinf(midline))
                 %        keyboard
                 %     end
                 %     midline = midline - mean(midline,1,"omitnan");
                 % 
                 %     % midline = lowpass(midline,4.0,data(i).Fs,'ImpulseResponse','fir');
                 %     % midline = lowpass(midline,0.01,data(i).Fs,'ImpulseResponse','fir');
                 % 
                 %     % midline = midline1.*midline2.*midline3.*midline4.*midline5.*midline6;
                 %     % midline = 0.5*midline/mean(midline);%normalize
                 % 
                 %     % midline = highpass(midline,.1,data(i).Fs,'ImpulseResponse','fir');
                 %     % midline_lp = lowpass(midline,.1,data(i).Fs,'ImpulseResponse','fir');
                 %     % midline_detrend = detrend_lin_pw(midline,357);
                 %     % 
                 %     % figure(4006), plot(midline), hold on, 
                 %     % plot(midline1)
                 %     % plot(midline_detrend)
                 %     % plot(mean(midline.'),'k'), plot(median(midline.'),'r');
                 %     % xlim([1000,5000])
                 %     values_aux(end+1) = values_aux(1); %duplicate an aux
                 %     keys_aux(end+1) = {'NIRmidline'}; %change its name
                 %     values_aux{end}.description = 'NIRmidline';
                 %     values_aux{end}.link = data(i).probe.link(lst,:);
                 %     values_aux{end}.time = data(i).time;
                 %     % values_aux{end}.Fs = data(i).Fs; %error: Fs dependent
                 %     % property and is set automatically to the correct value
                 %     % values_aux{end}.resampled = true;
                 %     values_aux{end}.record_start_delay_s = 0;
                 %     values_aux{end}.stimulus_delay = 0;
                 %     values_aux{end}.data = midline;
                 %     % figure(321312), clf, plot(data(i).time,midline)
                 %     % drawnow
                 % % else
                 % %     midline = zeros(size(data(1).data,1),1);
                 % %     disp([num2str(i),' midline has already been removed, midline regressor does not exist... using zeroes!'])
                 % % 
                 % %     values_aux(end+1) = values_aux(1); %duplicate an aux
                 % %     keys_aux(end+1) = {'NIRmidline'}; %change its name
                 % %     values_aux{end}.description = 'NIRmidline';
                 % %     % values_aux{end}.link = data(i).probe.link(lst,:);
                 % %     values_aux{end}.time = data(i).time;
                 % %     % values_aux{end}.Fs = data(i).Fs; %error: Fs dependent
                 % %     % property and is set automatically to the correct value
                 % %     % values_aux{end}.resampled = true;
                 % %     values_aux{end}.record_start_delay_s = 0;
                 % %     values_aux{end}.stimulus_delay = 0;
                 % %     values_aux{end}.data = midline;
                 % 
                 % end

                 %% keep only certain auxiliary data
                 for aa=1:numel(indexaux2keep)                     
                     indexaux2keep(aa) = find(contains(keys_aux,obj.aux2keep{aa},'IgnoreCase',true));
                 end
                 valuesaux_pruned = values_aux(indexaux2keep);
                 keysaux_pruned = keys_aux(indexaux2keep);

                 if obj.resampleaux %&& data(i).Fs~=auxFs
                     Fs_resample = (data(i).Fs);
                     multiplier = 1000; % workaround to be able to match the fNIRS sampling rate 8.93 Hz; 'resample' function requires integer ratios.
                     for vv=1:length(keysaux_pruned)
                         if valuesaux_pruned{vv}.Fs == Fs_resample
                             continue
                         else
                             Fs_i = round(valuesaux_pruned{vv}.Fs); %round() is necessary to avoid the following Error using resample: Expected input number 3, Q, to be integer-valued. Error in signal.internal.resample.validateResampleRatio  
                            % Fs_i = valuesaux_pruned{vv}.Fs;
                             tvec_resampled = valuesaux_pruned{vv}.time;
                             data_resampled = valuesaux_pruned{vv}.data;
                             if obj.trimtimevectors
                                 tmax = min(data(i).time(end),tvec_resampled(end));
                                 index_todiscard_attheend = nnz(tvec_resampled<=tmax)+1;
                                 data_resampled(index_todiscard_attheend:end)=[];
                                 % tvec_resampled(index_todiscard_attheend:end)=[];
                                 tmin = max(data(i).time(1),tvec_resampled(1));
                                 index_todiscard_atthebeginning = nnz(tvec_resampled<=tmin)-1;
                                 data_resampled(1:index_todiscard_atthebeginning)=[];
                                 % tvec_resampled(1:index_todiscard_atthebeginning)=[];
                             end

                             zz = data_resampled;
                             data_resampled = resample(data_resampled,round(Fs_resample*multiplier),Fs_i);
                             data_resampled = resample(data_resampled,1,multiplier);
                             while length(data_resampled)>size(data(i).data,1)
                                 if abs(length(data_resampled)-size(data(i).data,1))>10
                                    keyboard
                                 end
                                 if abs(tvec_resampled(end)-data(i).time(end)) < abs(tvec_resampled(1)-data(i).time(1))
                                     data_resampled(1)=[];
                                 else
                                     data_resampled(end) = [];
                                     % i
                                     % mfilename
                                 end
                             end
                             if length(data_resampled)~=size(data(i).data,1) %This should not happen if resampling is done properly
                                 keyboard
                             end
                             tvec_resampled = data(i).time; %make it identical to the fnirs time vector
                             % figur('resample'); clf; data(i).draw; hold on; plot(tvec_resampled,data_resampled)
                             valuesaux_pruned{vv}.data = data_resampled;
                             valuesaux_pruned{vv}.time = tvec_resampled;
                         end
                     end
                 end

                 if obj.visualize
                     figure(460000+i),clf,
                     N_aux = length(keysaux_pruned);
                     axpos=subplotpos(N_aux+1,1,'ygap',0);
                     ax=gobjects(N_aux+1,1);
                     for vv=1:N_aux
                         ax(vv)=axes;%subplot(N_aux+1,1,vv);
                         pl = plot(ax(vv),valuesaux_pruned{vv}.time,valuesaux_pruned{vv}.data);
                         pl.DisplayName = valuesaux_pruned{vv}.description;
                         hold on
                         legend('Interpreter','none')
                         ax(vv).Position= axpos(vv,:);
                         % plot(valuesaux_pruned{1}.time,valuesaux_pruned{1}.data),
                         % legend(valuesaux_pruned{vv}.description,valuesaux_pruned{1}.description)
                     end
                     ax(vv+1)=axes;
                     plot(ax(vv+1),data(i).time,data(i).data)
                     ax(vv+1).Position= axpos(vv+1,:);
                     linkaxes(ax,'x')
                     % legend('Interpreter','none')
                 end
                aux_pruned = Dictionary(keysaux_pruned,valuesaux_pruned);
                data(i).auxillary = aux_pruned;
            end
        end
    end
end