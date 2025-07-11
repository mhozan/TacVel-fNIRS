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
        aux2keep = {'ADIpulsetransducer',...
            'mic',...
            'Eyeblink_rate',...
            'EMG_face',...
            'hand',... %EMG2_hand or EMG2_righthand
            'XYZ_acceleration',...
            'Respiration_rate',...
            'CombinedAuxMotion',...
            'NIRmidline'};
        % aux2keep = {'XYZ_acceleration'};
        % resampleaux=true;
        resampleaux=true;
        trimtimevectors = true;
        clip_threshold = 99;

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

                % auxmultiplied = ones(size(values_aux(1)));


                %% accelerometer

                acclindices = find(contains(keys_aux,'accl','IgnoreCase',true));
                if ~isempty(acclindices)
                    combinedaccelerometer = values_aux(acclindices);
                    values_aux(end+1) = values_aux(acclindices(1)); %duplicate Xaccl
                    keys_aux(end+1) = {'XYZ_acceleration'}; %change its name
                    values_aux{end}.description = 'XYZ_acceleration';
                    
                    % thresh = .9;
                    lowerprctrile = 90;
                    multiplyfactor = 10;

                    thresh1low = prctile(combinedaccelerometer{1}.data,lowerprctrile);
                    thresh1hi = thresh1low*multiplyfactor;                   
                    combinedaccelerometer{1}.data(combinedaccelerometer{1}.data>thresh1hi)=thresh1hi;
                    combinedaccelerometer{1}.data(combinedaccelerometer{1}.data<thresh1low)=thresh1low;
                    combinedaccelerometer{1}.data = combinedaccelerometer{1}.data - thresh1low;
                    combinedaccelerometer{1}.data = combinedaccelerometer{1}.data./max(combinedaccelerometer{1}.data);
                    combinedaccelerometer{1}.data(isnan(combinedaccelerometer{1}.data)) = 1;

                    thresh2low = prctile(combinedaccelerometer{2}.data,lowerprctrile);
                    thresh2hi = thresh2low*multiplyfactor;                   
                    combinedaccelerometer{2}.data(combinedaccelerometer{2}.data>thresh2hi)=thresh2hi;
                    combinedaccelerometer{2}.data(combinedaccelerometer{2}.data<thresh2low)=thresh2low;
                    combinedaccelerometer{2}.data = combinedaccelerometer{2}.data - thresh2low;
                    combinedaccelerometer{2}.data = combinedaccelerometer{2}.data./max(combinedaccelerometer{2}.data);
                    combinedaccelerometer{2}.data(isnan(combinedaccelerometer{2}.data)) = 1;

                    thresh3low = prctile(combinedaccelerometer{3}.data,lowerprctrile);
                    thresh3hi = thresh3low*multiplyfactor;                   
                    combinedaccelerometer{3}.data(combinedaccelerometer{3}.data>thresh3hi)=thresh3hi;
                    combinedaccelerometer{3}.data(combinedaccelerometer{3}.data<thresh3low)=thresh3low;
                    combinedaccelerometer{3}.data = combinedaccelerometer{3}.data - thresh3low;
                    combinedaccelerometer{3}.data = combinedaccelerometer{3}.data./max(combinedaccelerometer{3}.data);                                       
                    combinedaccelerometer{3}.data(isnan(combinedaccelerometer{3}.data)) = 1;
                    
                    % thresh2 = prctile(combinedaccelerometer{2}.data,50);
                    % combinedaccelerometer{2}.data(combinedaccelerometer{2}.data>thresh2)=thresh2;
                    % thresh3 = prctile(combinedaccelerometer{3}.data,50);
                    % combinedaccelerometer{3}.data(combinedaccelerometer{3}.data>thresh3)=thresh3;
                    % combinedaccelerometerdata = [combinedaccelerometer{1}.data,...
                    %     combinedaccelerometer{2}.data,...
                    %     combinedaccelerometer{3}.data];

                    %multiply X Y Z axes of the accelerometer
                    % xyzaccl = combinedaccelerometer{1}.data.*...
                    %     combinedaccelerometer{2}.data.*...
                    %     combinedaccelerometer{3}.data;
                    % xyzaccl = abs(combinedaccelerometer{1}.data)+...
                    %     abs(combinedaccelerometer{2}.data)+...
                    %     abs(combinedaccelerometer{3}.data);   
                    xyzaccl = sqrt(... %Vector Magnitude aka Euclidian norm
                        combinedaccelerometer{1}.data.^2+...
                        combinedaccelerometer{2}.data.^2+...
                        combinedaccelerometer{3}.data.^2);   
                    xyzaccl = xyzaccl - median(xyzaccl);
                    % accl_normalizationfactor  = prctile(xyzaccl,99.9);
                    % xyzaccl = 1*xyzaccl/max(xyzaccl);
                    % xyzaccl = 1*xyzaccl/accl_normalizationfactor;

                    % Perform PCA
                    % [coeff, score, ~] = pca([...
                    %     combinedaccelerometer{1}.data,...
                    %     combinedaccelerometer{2}.data,...
                    %     combinedaccelerometer{3}.data]);
                    % % Extract the first principal component (dominant movement)
                    % xyzaccl = score(:,1);

                    % threshlow = prctile(xyzaccl,50);
                    % threshhigh = 10*threshlow;
                    % 
                    % 
                    % xyzaccl = highpass(xyzaccl,0.5,values_aux{acclindices(1)}.Fs,'ImpulseResponse','fir');
                    % xyzaccl = lowpass(xyzaccl,20,values_aux{acclindices(1)}.Fs,'ImpulseResponse','fir');
                    % xyzaccl = abs(xyzaccl);
                    % [yupper,~] = envelope(xyzaccl);


                    if any(isnan(xyzaccl)) || any(isinf(xyzaccl)) % || ~isnumeric(xyzaccl) || isempty(xyzaccl)
                        keyboard
                    end
                    winlen = round(0.1*combinedaccelerometer{1}.Fs);
                    % [yupper,~] = envelope(xyzaccl,winlen,'analytic');
                    [xyz_acc_env,~] = envelope(xyzaccl,winlen,'rms'); %rms envelope
                    % yupper = yupper - mean(yupper);
                    % yupper = yupper - min(yupper);
                    % % xyzaccl2= xyzaccl;
                    % % xyzaccl2((round(.4*end):round(.6*end)))=0;
                    % % xyzaccl2((1:round(.05*end)))=0;
                    % % xyzaccl2(round(.95*end):end)=0;
                    % % thresh2 = prctile(xyzaccl2,obj.clip_threshold);
                    % % thresh2 = accl_normalizationfactor;
                    % tobecappedindices = xyzaccl>threshhigh;
                    % xyzaccl(tobecappedindices) = threshhigh; %ceil the high values
                    % tobezeroedindices= xyzaccl<threshlow;
                    % % xyzaccl(xyzaccl<prctile(xyzaccl,threshlow)) = threshlow;
                    % xyzaccl(tobezeroedindices) = threshlow;
                    % % xyzaccl = 1*xyzaccl/thresh2;
                    % xyzaccl = xyzaccl/max(xyzaccl);
                    % xyzaccl = xyzaccl - min(xyzaccl);
                    % % xyzaccl = highpass(xyzaccl,5,values_aux{acclindices(1)}.Fs,'ImpulseResponse','fir');
                    % % xyzaccl = xyzaccl.^2;
                    % % xyzaccl(xyzaccl>prctile(xyzaccl,obj.clip_threshold)) = prctile(xyzaccl,obj.clip_threshold);
                    % % thresh = prctile(xyzaccl,99.99);
                    % % xyzaccl(xyzaccl>thresh)=thresh;
                    % % xyzaccl = 1*xyzaccl/thresh2;
                    % % xyzaccl = xyzaccl-mean(xyzaccl);
                    % % figure, plot(values_aux{acclindices(1)}.time,xyzaccl)


                    figur(['xyz_accl_',num2str(i,'%0.2d')]), clf
                    l1 = plot(combinedaccelerometer{1}.data); l1.DisplayName = 'x_acc';
                    hold on
                    l2 = plot(combinedaccelerometer{2}.data); l2.DisplayName = 'y_acc';
                    l3 = plot(combinedaccelerometer{3}.data); l3.DisplayName = 'z_acc';
                    l4 = plot(xyzaccl+1); l4.DisplayName = 'xyz_acc';
                    l5 = plot(xyz_acc_env+4,'r','LineWidth',2); l5.DisplayName = 'upper_env_analytic';
                    % l6 = plot(yupper2,'--k'); l6.DisplayName = 'upper_env_rms';
                    % ylim([0,0.5])
                    % figure, plot(xyzaccl)
                    legend('Interpreter','none')
% figure, plot(yupper)


                    % xyzaccl = xyzaccl/(max(xyzaccl));
                    % values_aux{end}.data = xyzaccl; %replace its data
                    values_aux{end}.data = xyz_acc_env; %replace its data
                    % auxmultiplied =  auxmultiplied.*xyzaccl;
                else
                    fprintf('no acceleromter auxiliary data found.')
                end

                %% Respiration
                respindices = find(contains(keys_aux,{'resp','cage','abdomen'},'IgnoreCase',true));
                if ~isempty(respindices)
                    RESP_aux = values_aux(respindices);
                    RESP_aux_combined = [RESP_aux{1}.data,...
                        RESP_aux{2}.data,RESP_aux{3}.data,...
                        RESP_aux{4}.data];
                    normalizeeachrespcolumn = true;
                    if normalizeeachrespcolumn
                        respnormalizationfactor  = prctile(RESP_aux_combined,90);
                        if all(respnormalizationfactor) % all resp channels are mostly non zero
                        RESP_aux_combined = RESP_aux_combined./respnormalizationfactor;
                        else
                        RESP_aux_combined = 0*RESP_aux_combined;
                        disp(['i = ',num2str(i,'%0.2d'),'Respiration rate data is not good. set to zero.'])
                        end
                    end

                    RESP_median = median(RESP_aux_combined,2)./60;

                    % figur(['respiration_',num2str(i,'%0.2d')]), clf
                    % plot(RESP_aux_combined), hold on
                    % plot(60*RESP_median)
                    % legend
                    % RESP_median = lowpass(RESP_median,20,values_aux{respindices(1)}.Fs,'ImpulseResponse','fir');
                    % RESP_median = highpass(RESP_median,.1,values_aux{respindices(1)}.Fs,'ImpulseResponse','fir');

                    if any(isnan(RESP_median))
                        keyboard
                    end

                    values_aux(end+1) = values_aux(respindices(1)); %duplicate RCLo_Resp
                    keys_aux(end+1) = {'Respiration_rate'}; %change its name
                    values_aux{end}.description = 'Respiration_rate';
                    values_aux{end}.data = RESP_median; %replace its data

                    % RESP_median(RESP_median>prctile(RESP_median,obj.clip_threshold)) = prctile(RESP_median,obj.clip_threshold);
                    % auxmultiplied =  auxmultiplied.*RESP_median/(max(RESP_median));



                else
                    fprintf('no resp auxiliary data found.')
                end

                %% Separate eyeblink and facial muscle tone
                eyeEMGindex = find(contains(keys_aux,'eye','IgnoreCase',true));
                if ~isempty(eyeEMGindex)
                    eyeEMGdata = values_aux{eyeEMGindex}.data;
                    eyeblink = eyeEMGdata(:,1);
                    facialEMG = (eyeEMGdata(:,2).^2);

                    values_aux(end+1) = values_aux(eyeEMGindex); %duplicate 'CH05_07_EMG1_eyelid'
                    keys_aux(end+1) = {'Eyeblink_rate'}; %change its name
                    values_aux{end}.description = 'Eyeblink_rate';
                    values_aux{end}.data = eyeblink; %replace its data

                    values_aux(end+1) = values_aux(eyeEMGindex); %duplicate 'CH05_07_EMG1_eyelid'
                    keys_aux(end+1) = {'EMG_face'}; %change its name
                    values_aux{end}.description = 'EMG_face';
                    facialEMG2= facialEMG;
                    facialEMG2((round(.4*end):round(.6*end)))=0;
                    facialEMG2((1:round(.05*end)))=0;
                    facialEMG2(round(.95*end):end)=0;
                    threshemg = prctile(facialEMG2,obj.clip_threshold);
                    facialEMG(facialEMG2>threshemg) = threshemg;
                    facialEMG = facialEMG./threshemg;
                    facialEMG = highpass(facialEMG,5,values_aux{eyeEMGindex}.Fs,'ImpulseResponse','fir');
                    facialEMG = facialEMG.^2;
                   
                    values_aux{end}.data = facialEMG; %replace its data
                    % auxmultiplied =  auxmultiplied.*facialEMG/(max(facialEMG));



                else
                    fprintf('no eye_emg auxiliary data found.')
                end

                %% CombinedAuxMotion: a multiplication of a few auxiliary channels to best capture motion artifacts in a single channel
                % clip_threshold = 99;
                 emghandindex = find(contains(keys_aux,'hand','IgnoreCase',true));
                 emghand = values_aux{emghandindex}.data;
                 emghand(emghand>prctile(emghand,obj.clip_threshold)) = prctile(emghand,obj.clip_threshold);
                 % auxmultiplied = auxmultiplied.*emghand/(max(emghand));

                 micindex = find(contains(keys_aux,'mic','IgnoreCase',true));
                 micdata = values_aux{micindex}.data;
                 micdata = highpass(micdata,10,values_aux{micindex}.Fs,'ImpulseResponse','fir');
                 micdata(micdata>prctile(micdata,obj.clip_threshold)) = prctile(micdata,obj.clip_threshold);
                 % auxmultiplied = auxmultiplied.*micdata/(max(micdata));

                 % thresh = prctile(auxmultiplied,obj.clip_threshold);
                 % auxmultiplied(auxmultiplied>thresh)=thresh;
                 % auxmultiplied=auxmultiplied./thresh;
                 fc_hipass= 1.5;
                 % auxmultiplied = highpass(auxmultiplied,fc_hipass,values_aux{1}.Fs,'ImpulseResponse','fir');
                 
 
                 values_aux(end+1) = values_aux(1); %duplicate an aux
                 keys_aux(end+1) = {'CombinedAuxMotion'}; %change its name
                 values_aux{end}.description = 'CombinedAuxMotion';
                 % values_aux{end}.data = auxmultiplied; %replace its data'
                 % combinedmotion =  sqrt(0.5*(facialEMG+xyzaccl).^2); %replace its data;
                 combinedmotion =  facialEMG+xyzaccl+micdata+emghand; %replace its data;
                 % combinedmotion = highpass(combinedmotion,.5,values_aux{1}.Fs,'ImpulseResponse','fir');
                 combinedmotion = combinedmotion.^2;
                 combinedmotion = highpass(combinedmotion,10,values_aux{1}.Fs,'ImpulseResponse','fir');
                 combinedmotion = abs(combinedmotion);
                 combinedmotion = combinedmotion - min(combinedmotion);
                 combinedmotion((combinedmotion>prctile(combinedmotion,90)))  = prctile(combinedmotion,90);

                 values_aux{end}.data = combinedmotion; %replace its data

                 %% Consider centered optode links (S8-D9, S9-D9, S9-D10) as regressors

                 lst = ismember(data(i).probe.link.source,[8,9])  & ismember(data(i).probe.link.detector,[9,10]);
                 if isempty(lst) || ~nnz(lst)
                     keyboard
                 end
                 if ~isempty(lst) && nnz(lst)
                     midline = data(i).data(:,lst);
                     midline1= midline(:,1);
                     midline2= midline(:,2);
                     midline3= midline(:,3);
                     midline4= midline(:,4);
                     midline5= midline(:,5);
                     midline6= midline(:,6);
                     % midline = mean(midline.').'; %unstable model?
                     midline = median(midline.').';
                     midline = midline-min(midline)+eps; %make sure deta is greater than zero

                     %optical density
                     m = mean( midline, 1 ,"omitnan");
                     midline =log(midline./(ones(size(midline,1),1)*m));
                     midline(midline<0)=0;
                     if any(isnan(midline)) || any(isinf(midline))
                        keyboard
                     end
                     midline = midline - mean(midline,1,"omitnan");

                     % midline = lowpass(midline,4.0,data(i).Fs,'ImpulseResponse','fir');
                     % midline = lowpass(midline,0.01,data(i).Fs,'ImpulseResponse','fir');

                     % midline = midline1.*midline2.*midline3.*midline4.*midline5.*midline6;
                     % midline = 0.5*midline/mean(midline);%normalize

                     % midline = highpass(midline,.1,data(i).Fs,'ImpulseResponse','fir');
                     % midline_lp = lowpass(midline,.1,data(i).Fs,'ImpulseResponse','fir');
                     % midline_detrend = detrend_lin_pw(midline,357);
                     % 
                     % figure(4006), plot(midline), hold on, 
                     % plot(midline1)
                     % plot(midline_detrend)
                     % plot(mean(midline.'),'k'), plot(median(midline.'),'r');
                     % xlim([1000,5000])
                     values_aux(end+1) = values_aux(1); %duplicate an aux
                     keys_aux(end+1) = {'NIRmidline'}; %change its name
                     values_aux{end}.description = 'NIRmidline';
                     values_aux{end}.link = data(i).probe.link(lst,:);
                     values_aux{end}.time = data(i).time;
                     % values_aux{end}.Fs = data(i).Fs; %error: Fs dependent
                     % property and is set automatically to the correct value
                     % values_aux{end}.resampled = true;
                     values_aux{end}.record_start_delay_s = 0;
                     values_aux{end}.stimulus_delay = 0;
                     values_aux{end}.data = midline;
                     % figure(321312), clf, plot(data(i).time,midline)
                     % drawnow
                 % else
                 %     midline = zeros(size(data(1).data,1),1);
                 %     disp([num2str(i),' midline has already been removed, midline regressor does not exist... using zeroes!'])
                 % 
                 %     values_aux(end+1) = values_aux(1); %duplicate an aux
                 %     keys_aux(end+1) = {'NIRmidline'}; %change its name
                 %     values_aux{end}.description = 'NIRmidline';
                 %     % values_aux{end}.link = data(i).probe.link(lst,:);
                 %     values_aux{end}.time = data(i).time;
                 %     % values_aux{end}.Fs = data(i).Fs; %error: Fs dependent
                 %     % property and is set automatically to the correct value
                 %     % values_aux{end}.resampled = true;
                 %     values_aux{end}.record_start_delay_s = 0;
                 %     values_aux{end}.stimulus_delay = 0;
                 %     values_aux{end}.data = midline;

                 end

                 %% keep only certain auxiliary data
                 for aa=1:numel(indexaux2keep)
                     indexaux2keep(aa) = find(contains(keys_aux,obj.aux2keep{aa},'IgnoreCase',true));
                 end
                 valuesaux_pruned = values_aux(indexaux2keep);
                 keysaux_pruned = keys_aux(indexaux2keep);

                 % auxFs = data(i).Fs;
                 % 
                 % if ~isempty(valuesaux_pruned)
                 %     auxFs = valuesaux_pruned{1}.Fs;
                 % end





                 if obj.resampleaux %&& data(i).Fs~=auxFs
                     Fs_resample = (data(i).Fs);
                     multplier = 1000; % workaround to be able to match the fNIRS sampling rate 8.93 Hz; 'resample' function requires integer ratios.
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
                                 index_todiscardfrom_attheend = nnz(tvec_resampled<tmax);
                                 % tvec_resampled(index_todiscardfrom_attheend+1:end)=[];
                                 data_resampled(index_todiscardfrom_attheend+1:end)=[];
                             end
                             % tvec_resampled = resample(tvec_resampled,round(Fs_resample*multplier),Fs_i);
                             % tvec_resampled = resample(tvec_resampled,1,multplier);
                             % tvec_resampled = valuesaux{vv}.time(1)+raw.time;

                             % tvec_resampled(tvec_resampled>(max(raw.time)))=[];
                             zz = data_resampled;
                             data_resampled = resample(data_resampled,round(Fs_resample*multplier),Fs_i);
                             data_resampled = resample(data_resampled,1,multplier);
                             % data_resampled = resample(data_resampled,multplier,round(Fs_resample*multplier));

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
                             if length(data_resampled)~=size(data(i).data,1)
                                 keyboard
                             end

                             tvec_resampled = data(i).time; %make it identical to the fnirs time vector


                             % figur('resample'); clf; data(i).draw; hold on; plot(tvec_resampled,data_resampled)
                             valuesaux_pruned{vv}.data = data_resampled;
                             valuesaux_pruned{vv}.time = tvec_resampled;
                         end
                         % valuesaux{vv}.Fs = Fs_resample;%automatically will be set

                     end

                 end

                 figure(460000+i),clf,
                 for vv=1:length(keysaux_pruned)
                     pl = plot(valuesaux_pruned{vv}.time,valuesaux_pruned{vv}.data);
                     pl.DisplayName = valuesaux_pruned{vv}.description;

                     hold on
                     % plot(valuesaux_pruned{1}.time,valuesaux_pruned{1}.data),
                     % legend(valuesaux_pruned{vv}.description,valuesaux_pruned{1}.description)
                 end
                 legend('Interpreter','none')
% if i==21
%     keyboard
% end
                aux_pruned = Dictionary(keysaux_pruned,valuesaux_pruned);
                data(i).auxillary = aux_pruned;

            end
        end
    end
end









