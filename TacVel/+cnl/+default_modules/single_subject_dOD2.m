function jobs = single_subject_dOD2
jobs=nirs.modules.ImportData();
jobs.Input='raw';
jobs=nirs.modules.RemoveStimless(jobs);

jobs = nirs.modules.FixNaNs(jobs);
jobs = nirs.modules.Resample(jobs);
jobs.Fs = 2; % resample to 5 Hz

jobs = nirs.modules.OpticalDensity( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='dOD';
% jobs = nirs.modules.TrimBaseline( jobs );
% jobs.preBaseline   = 30;
% jobs.postBaseline  = 30;

jobs = nirs.modules.BeerLambertLaw( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='Hb';

jobs = nirs.modules.GLM(jobs );
jobs.type='AR-IRLS'; %default
jobs.AddShortSepRegressors=true;
jobs = nirs.modules.AddAuxRegressors(jobs);

jobs = nirs.modules.ExportData(jobs);
jobs.Output='SubjStats';

jobs = nirs.modules.BaselineCorrection(jobs);
jobs.Output='Baselinecorrected';




