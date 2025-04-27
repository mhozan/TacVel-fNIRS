function [] = expfit(f_vec,spctrgrm)%[expcurv,expcurv_inv] = expfit(f_vec,spctrgrm)
% [expcurv,expcurv_inv] = EXPFIT(f_vec,spctrgrm); expcurv is a vector the same size
% as f_vec which is a two-term exponential fit to the intensity profile of the
% spctrgrm. spctgrm is an MxN matrix where M must be equal to the length of
% f_vec. 
%      General model Exp2:
%      f2(x) =  ;
%mohsen hozan Oct 2016

intensity_profile = sum(spctrgrm');
intensity_profile = intensity_profile./max(intensity_profile);
[~,peak_index] = max(intensity_profile);
f_peak=f_vec(peak_index);


x = f_vec(:);
y = intensity_profile(:);

% exponfit = fittype( @(a, b, c, d, e, x) a*exp(b*x) + c*exp(d*x) + e , 'problem' , {'a', 'b', 'c', 'd', 'e'}, 'independent', 'x', 'dependent', 'z')
exponfit = fittype( @(a, b, f01, c, d, f02, e, x) a*exp(b*(x+f01)) + c*exp(d*(x+f02)) + e)



%%
% f = fit(x,y,exponfit)
f = fit(x,y,exponfit,'Startpoint',[1 -1 f_peak -.1 .1 0 0])
x2=-1:0.1:30;
y_fit = f.a*exp(f.b*(x2+f.f01)) + f.c*exp(f.d*(x2+f.f02)) + f.e;
figure(678746), clf, 
plot(x,y,'.-.r') %the intensity curve
hold on
plot(x2,y_fit)
legend('data','exp2 fit')


% at

% expvec = exp(-alfa*x_vec())







