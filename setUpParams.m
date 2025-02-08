function params = setUpParams(model)

% Function to set defaults for recording and regression parameters 

% model: which model to run
%       model1: b-spline basis set kernels, fit each animal and day separately; option for regularization 
%       model2: fit (for mixed-effects, ME=1, only if regFlag == 0) regression concatenating all animals; b-spline basis set kernels; fit each session separately; option for regularization
%       model3: fit regression by animal, concatenating across days; b-spline basis set kernels, option for regularization and mixed effects (only if regFlag == 0) 

% Data location
if ispc
    params.fbasename = 'C:\Users\jmc0163\OneDrive - Northwestern University\Data\DA_Photometry_GLopez';
    params.fbasename_raw = 'C:\Users\jmc0163\OneDrive - Northwestern University\GLopez_DataforRegressionAnalysis';
elseif ismac
    params.fbasename = '/Users/julia/Library/CloudStorage/OneDrive-NorthwesternUniversity/Data/DA_Photometry_GLopez';
    params.fbasename_raw = '/Users/julia/Library/CloudStorage/OneDrive-NorthwesternUniversity/GLopez_DataforRegressionAnalysis';
else
    params.fbasename = '/projects/p31438/DA_Photometry_GLopez/';
end

% Regression parameters 
params.model        = model;
switch model
     case 'model1' % fit by animal and day; basis set
        % Fit parameters
        params.eventNames   = {'CueAvoid';'CueEscape';'Shock';'AvoidCross';'EscapeCross'}; % which events to include
        params.timeBack     = [2;       2;          .5;        1;            0]; % time before the event in seconds
        params.timeForward  = [1.5;       1.5;        1.5;        5;            5]; % time after the event in seconds
        params.sessIDs      = 1:7; % which sessions to fit
        params.speedFlag    = false;
    case 'model2' % fit by session; concatenate all mice
        % Fit parameters
        params.eventNames   = {'CueAvoid';'CueEscape';'Shock';'AvoidCross';'EscapeCross'}; % which events to include
        params.timeBack     = [2;       2;          .5;        1;            0]; % time before the event in seconds
        params.timeForward  = [1.5;      1.5;        1.5;        5;            5]; % time after the event in seconds
        params.sessIDs      = 1:7; % which sessions to fit
        params.speedFlag    = false;
    case 'model3' % includes speed 
        % Fit parameters 
        params.eventNames   = {'CueAvoid';'CueEscape';'Shock';'AvoidCross';'EscapeCross';'ITICross'}; % which events to include
        params.timeBack     = [2;       2;          .5;        1;            0;           2]; % time before the event in seconds
        params.timeForward  = [1.5;      1.5;        1.5;        5;          5;           5]; % time after the event in seconds
        params.sessIDs      = 1:7; % which sessions to fit
        params.speedFlag    = true; % include speed predictor 
end

params.numBasis     = 7; % number of functions per second for basis set
params.nFolds       = 5; % number of folds for cross validation; if == 1, uses leave one out cross validation
params.crossValFlag = 1; % perfom cross validation? 
params.regFlag      = 0; % perform regularization? 0 = no, 1 = ridge, 2 = lasso
params.MEFlag       = 0;

% Shuffles 
params.numShuff = 0; %number of resamples for bootstrapping

% Data parameters
params.newFs          = 20; % frequency to downsample to
params.Fs             = round(1017.25); % acquisition rate
params.eventNames_all = {'Shock';'EscapeCross';'AvoidCross';'CueEscape';'CueAvoid';'ITICross'}; % all events
params.regions        = {'vmShell';'Core'}; % brain regions
params.numBasis       = params.numBasis;

% Plotting colors
params.avoidC     =[32 118 188]./255;
params.escapeC    =[242 103 25]./255;
params.shockC     = [127 22 26]./255;

% Which kernels to compare 
params.whichTests = false(numel(params.eventNames),numel(params.eventNames));
params.whichTests((params.eventNames == "CueAvoid"),(params.eventNames=="CueEscape")) = true; % compare cue avoid and cue escape kernels

%% Find median time between events

params = getLatencies_activeAvoid(params,params.fbasename);

for nr = 1:numel(params.regions)
    eval(sprintf('params.escapeMed.%s = median(cell2mat(params.escapeLatency.%s''));',params.regions{nr},params.regions{nr}));
    eval(sprintf('params.avoidMed.%s = median(cell2mat(params.avoidLatency.%s''));',params.regions{nr},params.regions{nr}));
end
