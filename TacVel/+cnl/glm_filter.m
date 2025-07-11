%glmval


% raw_all(7).draw
%  Subj_Stats = j.run(Hb(7));
% hrf = Subj_Stats.HRF
% 
% 
function yhat = glm_filter(b_coefficients, predictor_timeseries, auxiliary_timeseries, HRFs)
    % Input:
    %   b_coefficients: Coefficients of the predictors (including the intercept) [1 x num_predictors]
    %   predictor_timeseries: Predictor time series [num_samples x num_predictors]
    %   auxiliary_timeseries: Auxiliary time series to be regressed out [num_samples x num_auxiliary]
    %   HRFs: Canonical Hemodynamic Response Functions (HRFs) [hrf_length x num_predictors]
    %
    % Output:
    %   yhat: Filtered output time series (predicted signal without the influence of auxiliary timeseries) [num_samples x 1]

    % Ensure column vectors
    b_coefficients = b_coefficients(:);
    num_predictors = length(b_coefficients) - 1; % Subtract 1 for the intercept

    % Check input sizes
    [num_samples, num_auxiliary] = size(auxiliary_timeseries);
    assert(size(predictor_timeseries, 1) == num_samples, 'Predictor and auxiliary time series must have the same number of samples.');
    assert(size(HRFs, 1) == num_samples, 'Number of rows in HRFs must be equal to the number of samples.');

    % Create the design matrix with predictors and auxiliary timeseries
    X = [ones(num_samples, 1), predictor_timeseries, auxiliary_timeseries];

    % Apply regression to obtain the filtered output (predicted signal)
    yhat = X * b_coefficients;

    % Remove the contribution of auxiliary timeseries using GLM
    for i = 1:num_auxiliary
        yhat = yhat - conv(auxiliary_timeseries(:, i), HRFs(:, i), 'full');
    end

    % Only keep the original samples (discard the extra samples introduced by convolution)
    yhat = yhat(1:num_samples);
end
