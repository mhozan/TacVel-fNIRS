function [tvec_common,data1_synced,data2_synced]=synchronize2(tvec1,data1,tvec2,data2,f_common)
% uses matlab timeseries/synchronize to synchronize two time vectors and datasets to a new frequency, f_common (Hz)
% tvec1 is a 1xN vector; data1 is a RxN matrix. the sampling frequency
% across all rows of data1 must be consistent with tvec1.
% tvec2 is a 1xM vector; data2 is a PxM matrix
% EXAMPLE 
% % % % % % tvec1 = -.3:.1:1;
% % % % % % data1 = [tvec1.^2;tvec1-1];
% % % % % % tvec2 = .1:.13:1.5;
% % % % % % data2 = [0.5*tvec2+1;tvec2.^3];
% % % % % % f_common = 100; %(Hz) 0.01 interval
% % % % % % [tvec_common,data1_synced,data2_synced]=synchronize2(tvec1,data1,tvec2,data2,f_common);
% % % % % % figure(11); clf;
% % % % % % ax(1)=subplot(4,1,1); plot(tvec1,data1,'.'), hold on, title(ax(1),'data1')
% % % % % % ax(2)=subplot(4,1,2); plot(tvec2,data2,'.'), title(ax(2),'data2')
% % % % % % ax(3)=subplot(4,1,3); plot(tvec_common,data1_synced,'.'), title(ax(3),'data1:synced')
% % % % % % ax(4)=subplot(4,1,4); plot(tvec_common,data2_synced,'.'), title(ax(4),'data2:synced')
% % % % % % linkaxes(ax,'x'); xlim([min([tvec1(:);tvec2(:)]),   max([tvec1(:);tvec2(:)])])
% Mohsen 10/17/2016
% see alo SYNCHRONIZE, TIMESERIES CLASS
validateattributes(f_common,    {'numeric'},    {'scalar','nonempty','>',eps});
validateattributes(tvec1,       {'numeric'},    {'vector','nonempty'});
validateattributes(data1,       {'numeric'},    {'2d','nonempty'});
validateattributes(tvec2,       {'numeric'},    {'vector','nonempty'});
validateattributes(data2,       {'numeric'},    {'2d','nonempty'});
intervl = 1/f_common;
tol = 1e-10;
% data1_synced=[];
% data2_synced=[];

R = size(data1,1);
P = size(data2,1);

for r=1:R
    for p=1:P
        ts1=timeseries(data1(r,:),tvec1);
        ts2=timeseries(data2(p,:),tvec2);
        [ts1a,ts2a] = synchronize(ts1,ts2,'Uniform','Interval',intervl,'tolerance',tol);
        data2_synced(p,:) = ts2a.data;
%         disp(r)
%         disp(p)
    end
    data1_synced(r,:) = ts1a.data;
end
tvec_common = ts1a.Time;

