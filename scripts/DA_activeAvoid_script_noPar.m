%% See setUpParams.m for defaults



%% Set parameters
%%% Basis sets, fit by day, concatenated by animal, no regularization or cross validation, a few shuffles for testing 
model       = 'model1';
params              = setUpParams(model);
params.numShuff     = 0; 
params.MEFlag       = false;
params.regFlag      = 2;

plotFlag    = true; % plot?
fitFlag     = false; % fit model? 
extractData = false; % reorganize data? 
parFlag     = false; % use parfor?

%% Fit model 

if fitFlag 
    DA_activeAvoid(params,extractData,parFlag);
end

if plotFlag 
    plotModelFig(params)
end

