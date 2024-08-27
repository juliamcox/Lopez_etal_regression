function stats = plotModelFig(params)

%% Load data     
load(fullfile(params.fbasename,'plots',sprintf('%s_kernels_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag))),'temporalKernels','shuffKernels');
for nr = 1:numel(params.regions)
    saveName = fullfile(params.fbasename,'plots',params.regions{nr},sprintf('%s_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag)));
    load(saveName,'modelStats_real','modelStats_shuff')
    realStats{nr} = modelStats_real;
    shuffStats{nr} = modelStats_shuff;
end

%% Generate filename for saving figures 
fname = sprintf('modelFig_%s_Fs%s_numShuff%s_ME%s_reg%s',params.model,num2str(params.newFs),num2str(params.numShuff),num2str(params.MEFlag),num2str(params.regFlag));
saveLoc = fullfile(params.fbasename,'plots',params.model);
if ~isfolder(saveLoc)
    mkdir(saveLoc)
end

%% Plot temporal kernels (A-B)
CIFlag = true;
sessIDs = [1,3,7]; % which sessions to plot kernels
plotKernels(temporalKernels,sessIDs,params.timeBack,params.timeForward,params.eventNames,params,saveLoc,fname,shuffKernels,CIFlag);

%% Plot correlation coefficients for full model and comparison between full and events in eventNames (C-D)
if params.numShuff > 0
    CIFlag = true;
    stats = plotBootCorr(realStats,shuffStats,params,{'Shock'},CIFlag,saveLoc,fname);
end

%% Plot area under the curve (E)
sessIDs = 1:7;
CIFlag  = true;
kernelAUC(temporalKernels,shuffKernels,params,params.whichTests,{'Core'},sessIDs,CIFlag,saveLoc,fname);

