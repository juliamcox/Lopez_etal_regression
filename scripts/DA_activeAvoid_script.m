%% See setUpParams.m for defaults

% bash submit_MATLAB_job_regression.sh to submit job on Quest

%% Set parameters
%%% Basis sets, fit by day, concatenated by animal
% model       = 'model2';
% params              = setUpParams(model);
% params.numShuff     = 5000; 
% params.MEFlag       = false; 
% params.regFlag      = 2;
% 
% plotFlag    = true; % plot?
% fitFlag     = false; % fit model? 
% extractData = false; % reorganize data? 
% parFlag     = false; % use parfor for bootstrapping model2?
% 
% %% Set parameters
% %%% Basis sets, fit by day, and animal, no regularization or cross
% %%% validation, a few shuffles for testing ??????????
% model       = 'model1';
% params              = setUpParams(model);
% params.crossValFlag = 0; 
% params.MEFlag       = false; 
% params.regFlag      = 2;
% 
% plotFlag    = false; % plot?
% fitFlag     = true; % fit model? 
% extractData = false; % reorganize data? 
% parFlag     = true; % use parfor?

%% Set parameters
%%% Basis sets for events, continuous speed variable, fit by day, concatenated by animal, lasso regularization 
% model       = 'model3';
% params              = setUpParams(model);
% params.numShuff     = 5000; 
% params.MEFlag       = false; 
% params.regFlag      = 2;
% 
% plotFlag    = true; % plot?
% fitFlag     = false; % fit model? 
% extractData = false; % reorganize data? 
% parFlag     = true; % use parfor for bootstrapping model2?

%% Set parameters
%%% Basis sets for events, continuous speed and acceleration variables, fit by day, concatenated by animal, lasso regularization 
% model       = 'model4';
% params              = setUpParams(model);
% params.numShuff     = 20; 
% params.MEFlag       = false; 
% params.regFlag      = 2;
% 
% plotFlag    = true; % plot?
% fitFlag     = true; % fit model? 
% extractData = false; % reorganize data? 
% parFlag     = true; % use parfor for bootstrapping model2?

%% Set parameters
%%% Basis sets for events, continuous speed and acceleration variables, fit
%%% by day, concatenated by animal, lasso regularization; all data points
%%% (not just during events)
model       = 'model5';
params              = setUpParams(model);
params.numShuff     = 50; 
params.MEFlag       = false; 
params.regFlag      = 2;

plotFlag    = false; % plot?
fitFlag     = true; % fit model? 
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

