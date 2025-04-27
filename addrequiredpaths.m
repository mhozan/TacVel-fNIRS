function [] = addrequiredpaths()
%ADDREQUIREDPATHS adding the paths to folders such as chronux, etc.
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% March 2017
[FDpath,~,~]  = fileparts(mfilename('fullpath'));
C = strsplit(FDpath,'ForceDataClass'); oneup = C{1}; %one folder up to the Matlab Libraries

FDdefaultspath	= fullfile(FDpath,'FDdefaults');
FDtoolspath   	= fullfile(FDpath,'FDtools'); %,filesep
fnirsCNL        = fullfile(oneup,'fnirsCNL'); %,filesep
% chronuxpath     = genpath(fullfile(oneup,'Chronux_spectral_analysis')); %chronux 2.11--spectral analysis for multi taper spectrograms
utilgeneral     = genpath(fullfile(oneup,'utils_general',filesep)); %general utility m files such as printPDF,stringSEMICOLdelimitter, etc.  mainly written by myself, some from mathworks fileexhange
fdasrvfpath     = genpath(fullfile(oneup,'fdasrvf_MATLAB-master',filesep)); %fdasrvf_MATLAB-master m files for curve registration/alignment. written by Tucker et al.
% nirstoolbox     = genpath(fullfile(oneup,'huppertt-nirs-toolbox-d4cdc6c4c618',filesep)); %nirs-toolbox by Huppert lab.
% nirstoolbox     = genpath(fullfile(oneup,'huppertt-nirs-toolbox-20230418',filesep)); %nirs-toolbox by Huppert lab 2023 update.
nirstoolbox     = genpath(fullfile(oneup,'nirs-toolbox-master-20250416',filesep)); %nirs-toolbox by Huppert lab 2025 update.
edfread         = fullfile(oneup,'edfreadZip',filesep); %edfread: mathworks.com/matlabcentral/fileexchange/31900-edfread
adisdk          = fullfile(oneup,'JimHokanson-adinstruments_sdk_matlab-1061f17',filesep);
PACE_fda        = genpath(fullfile(oneup,'PACE_matlab-master','release2.17','PACE')); %http://www.stat.ucdavis.edu/PACE/
Homer3          = genpath(fullfile(oneup,'Homer3_master',filesep));
avolu_tCCA_GLM  = genpath(fullfile(oneup,'avolu_tCCA_GLM_master',filesep));
avolu_GLM_BCI   = genpath(fullfile(oneup,'avolu_GLM_BCI_master',filesep));
addpath(FDpath...
    ,FDdefaultspath...
    ,FDtoolspath...
    ,fnirsCNL...,chronuxpath...
    ,utilgeneral...
    ,fdasrvfpath...
    ,nirstoolbox...
    ,edfread...
    ,adisdk...
    ,PACE_fda...
    ,avolu_tCCA_GLM...
    ,avolu_GLM_BCI...
    ) ;


%restoredefaultpath


