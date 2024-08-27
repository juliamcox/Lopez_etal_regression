function stats = plotByEvent(params)
% Plot figure S3 for individual animal fit

%% Load data     
load(fullfile(params.fbasename,'plots',sprintf('%s_kernels_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag))),'temporalKernels');
for nr = 1:numel(params.regions)
    saveName = fullfile(params.fbasename,'plots',params.regions{nr},sprintf('%s_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag)));
    load(saveName,'modelStats_real')
    realStats{nr} = modelStats_real;
end

%% Generate filename for saving figures
fname = sprintf('modelFig_%s_Fs%s_numShuff%s_ME%s_reg%s',params.model,num2str(params.newFs),num2str(params.numShuff),num2str(params.MEFlag),num2str(params.regFlag));
saveLoc = fullfile(params.fbasename,'plots',params.model);
if ~isfolder(saveLoc)
    mkdir(saveLoc)
end

stats = multCompareCorr(realStats,params,saveLoc,fname);


