function binary_output = binarythresholder(sig,threshold_amp,threshold_dur)
% binaryoutput=BINARYTHRESHOLDER(sig,threshold_amp,threshold_dur)
% inputs:
%   sig:            input vector(normaly a moving_averaged output from NLX3.EESD or any other time series).
%                   if sig is already in binary format, use a threshold_amp of 0.5
%   theshold_amp:	the cutoff amplitude for the binary decision.  If sig is
%                   normalized to max value, a threshold like 0.1 is reasonable.
%   theshold_dur:	minimum duration(in samples) for a decision (0 or 1)to be considered
%                   valid, i.e. to ignore very sudden threshold passing. For instance, if sig has been above threshold_amp for a duration
%                   longer than threshold_dur, and then drops below threshold_amp for a few
%                   samples and then immediately goes back above threshold_amp,
%                   "threshold_dur" helps the function to ignore the duration below the threshold_amp.
% outputs:
%   binary_output:	a binary vector the same size as the input "sig"
%
%example:
% binaryoutput = binarythresholder([0 1 1 1 0 1 0 0 1 0 1 1 1 0 1 1 0 0 1 1 0 0 0 0 1],0.5,3);
% binaryoutput =                   [0 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]
%by Mohsen hozan@mit.edu 5/18/2016

initial_decision = single(sig(:)>threshold_amp); %conv does not support logical inputs x-(
conv_win = ones(threshold_dur,1); %to continousely count the number of ones(and zeros) in a window with length threshold_dur

sameeventcounts = conv(initial_decision,conv_win,'same');
% figure, plot(sameeventcounts)
% hold on , stairs([0 45000],[threshold_dur/2 threshold_dur/2],'r')


% z= medfilt1(z,500,'truncate');
% figure, plot(z)
% hold on , plot([0 45000],[threshold_dur/2 threshold_dur/2],'r')

% thresholdpassings = z>threshold_dur/2;
% count=0;
% for i=1:length(thresholdpassings)
%     if thresholdpassings(i)==1
%         count=count+1;
%     end
% end

intermediate_decision = sameeventcounts>threshold_dur/2;
temp3=0;
count=0;
% figure(3215), clf, ah(1)=subplot(2,1,1); stairs(intermediate_decision)
while ~isempty(temp3) && count<20  %getting rid of short_duration changes (the goal is to not find a shortdurationchange anymore(empty temp3). counter is only to prevent and infinity loop)
    temp0 = diff(intermediate_decision);
    temp1 = find(temp0);
    temp2 = diff(temp1);
    temp3 = find(temp2<threshold_dur);
    for i=1:length(temp3)
        indices_of_shortdurationchange = temp1(temp3(i))+1:temp1(temp3(i)+1);
        intermediate_decision(indices_of_shortdurationchange) = ~intermediate_decision(indices_of_shortdurationchange); %simply flip the shortdurationchange
    end
    count=count+1;
end
% figure(3215), ah(2)=subplot(2,1,2); stairs(intermediate_decision)
% linkaxes(ah)
% ylim([-0.1 1.1])

% for

binary_output = intermediate_decision;
% figure, stairs(binary_output)








