% See setUpParams.m for defaults

% bash submit_MATLAB_job_regression.sh to submit job on Quest

%% Set parameters
%%%%%%%% Model 2: Basis sets, fit by day, concatenated by animal %%%%%%%%
% model       = 'model2';
% params              = setUpParams(model);
% params.numShuff     = 5000; 
% params.MEFlag       = false; 
% params.regFlag      = 2; % lasso regularization
% 
% plotFlag    = true; % plot?
% fitFlag     = false; % fit model? 
% extractData = false; % reorganize data? 
% parFlag     = false; % use parfor for bootstrapping model2?
% 

%%%%%%%% Model 1: Basis sets, fit by day, and animal, no regularization or cross validation %%%%%%%%
model               = 'model1';
params              = setUpParams(model);
params.crossValFlag = 0; 
params.MEFlag       = false; 
params.regFlag      = 2;

plotFlag    = false; % plot?
fitFlag     = true; % fit model? 
extractData = false; % reorganize data? 
parFlag     = true; % use parfor?

%%%%%%%% Model 3: Basis sets for events, continuous speed variable, ITI cross event, fit by day, concatenated by animal, lasso regularization 
model               = 'model3';
params              = setUpParams(model);
params.numShuff     = 5000; 
params.MEFlag       = false; 
params.regFlag      = 2;

plotFlag    = true; % plot?
fitFlag     = false; % fit model? 
extractData = false; % reorganize data? 
parFlag     = true; % use parfor for bootstrapping model2?



%% Fit model 

if fitFlag 
    DA_activeAvoid(params,extractData,parFlag);
end

%% Plot 
if plotFlag 
    stats = plotModelFig(params);
end

