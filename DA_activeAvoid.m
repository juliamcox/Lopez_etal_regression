function DA_activeAvoid(params,extractData,parFlag)

% Fits multiple linear regression to photometry data using a b-spline basis set
% Saves goodness of fit measures, event triggered fluorescence and estimated temporal kernels for each event

% Inputs:
% params: data structure to set parameters for cross validation, shuffles, mixed effects, etc. See setUpParams.m for fields
% extractData: does the data need to be reformatted? (true or false)
% parFlag:

rng('shuffle'); % set random number generator seed based on current time

% Find animal IDs
recs = readtable(fullfile(params.fbasename,'Photometry_Identification.xlsx')); % Load list of recordings
recs.Properties.VariableNames = {'ID';'Fname';'Day';'Region'};
Fnames = recs.Fname;
Fnames = cellfun(@(x) x(2:end-1), Fnames, 'UniformOutput', false);
recs.Fname = string(Fnames);

% Excluded mice
excl = readtable(fullfile(params.fbasename,'ExcludedMice_MouseID.csv'),ReadVariableNames=false);
%excl.Var1(end+1) = {'152-071'};
%excl.Var1(end+1) = {'310-907'};

%% Reorganize/load data

if extractData
    if ~ispc
        error('Raw data not copied to Quest')
    else
        extractEscapeAvoid(params.fbasename_raw,params.fbasename,params);
    end
end


%% Fit model

switch params.model

    %% Model 1: Fit each animal and session separately; basis sets
    case "model1"
        for nr = 1:numel(params.regions)
            thisIDs = unique(recs.ID(contains(recs.Region,params.regions{nr})));
            thisIDs(contains(thisIDs,excl.Var1)) = [];
            for ns = params.sessIDs % fit data by session
                lastwarn(''); % reset warning indicator
                for na = 1:numel(thisIDs)

                    % Fit each animal separately
                    fprintf('Fitting %s mouse %s day %s\n', params.regions{nr}, thisIDs{na}, num2str(ns));
                    dataLoc = fullfile(params.fbasename,params.regions{nr},thisIDs{na}); % location of data for this mouse
                    try % Catch statement for 310-907 who doesn't have data on session 4

                        %%% Create predictor matrix and response vector
                        [X,Y,trialIndicator,numFun,basisSets,thisFluor,allFluor] = createPredictorMatrix_activeAvoid(params,params.eventNames,params.timeBack,params.timeForward,ns,dataLoc,0);
                        %%% Extract avg fluorescence for each event
                        for ne = 1:numel(params.eventNames)
                            eval(sprintf('meanFluor.%s.%s(na,:,ns) = thisFluor{ne};',params.regions{nr},params.eventNames{ne}));
                        end
                        %%% If fitting shuffles for null distribution, generate Y
                        for nss = 1:params.numShuff
                            [~,Y_shuff{nss},~,~,~,~,~] = createPredictorMatrix_activeAvoid(params,params.eventNames,params.timeBack,params.timeForward,ns,dataLoc,1);
                        end
                        %%% Fit full and reduced models and to real data and calculate F statistic for compariosn between full and all reduced models(reduced models exclude all predictors associated w/ a particular event)
                        if params.regFlag > 0
                            fprintf('Finding lambda value \n')
                            lam=calcLambda(params.regFlag,X,Y,ones(size(Y)),trialIndicator);
                        else
                            lam = NaN;
                        end

                        %%%% Fit each animal separately without cross validation
                        fprintf('Fitting regression \n')
                        [fits_real{na,ns},modelStats_real{na,ns},errFlag] = fitModel_activeAvoid(X,Y,params.eventNames,numFun,params.regFlag,ones(size(Y)),0,lam);

                        if errFlag
                            % If fitting throws a warning, exclude fit (everything to NaN)
                            modelStats_real{na,ns}.corr.full = NaN;
                            modelStats_real{na,ns}.corr_noRefit.full = NaN;
                        else
                            % Calculate correlation coefficient between real and estimated data for full and reduced models
                            if params.regFlag == 0
                                thisX = cat(2,ones(size(X,1),1),nanzscore(X,1)); % add column of ones to predictor matrix (intercept)
                                modelStats_real{na,ns} = calcModel_corr(thisX,nanzscore(Y,1),fits_real{na,ns},numFun, modelStats_real{na,ns},params.eventNames,params.regFlag,params.MEFlag);
                            else
                                modelStats_real{na,ns} = calcModel_corr(nanzscore(X,1),nanzscore(Y,1),fits_real{na,ns},numFun, modelStats_real{na,ns},params.eventNames,params.regFlag,params.MEFlag);
                            end
                        end

                        %%% Fit full and reduced models to shuffled data to generate null distribution of F-statisitc while accounting for temporal structure of the photometry data

                        if params.numShuff > 0
                            if errFlag % only fit if the fit to the real data worked
                                modelStats_real{na,ns}.pvals = NaN;
                            else
                                fprintf('Fitting shuffles');
                                for nss = 1:params.numShuff
                                    fprintf('.')
                                    [fits_shuff{na,ns}{nss},modelStats_shuff{na,ns}{nss},shuffErr]  = fitModel_activeAvoid(X,Y_shuff{nss},params.eventNames,numFun,params.regFlag,ones(size(Y)),0,lam);
                                    if ~shuffErr
                                        if params.regFlag == 0
                                            thisX = cat(2,ones(size(X,1),1),X);
                                            modelStats_shuff{na,ns}{nss} = calcModel_corr(thisX,Y,fits_shuff{na,ns}{nss},numFun, modelStats_shuff{na,ns}{nss},params.eventNames,params.regFlag,params.MEFlag);
                                        else
                                            modelStats_shuff{na,ns}{nss} = calcModel_corr(X,Y,fits_shuff{na,ns}{nss},numFun, modelStats_shuff{na,ns}{nss},params.eventNames,params.regFlag,params.MEFlag);
                                        end
                                    else
                                        modelStats_shuff{na,ns}{nss}.corr.full = NaN;
                                    end
                                end
                                fprintf('\n')
                                % Calculate p-values
                                if isnan(fits_real{na,ns}.betas.full) % if the real data didn't fit properly, don't calculate p-values
                                    modelStats_real{na,ns}.pvals = NaN;
                                else
                                    modelStats_real{na,ns} = getPvals(modelStats_real{na,ns},modelStats_shuff{na,ns},params.eventNames);
                                end
                            end
                        else
                            modelStats_real{na,ns}.pvals = NaN;
                            fits_shuff{na,ns} = NaN;
                            modelStats_shuff{na,ns} = NaN;
                        end


                        %%% Organize data for plotting
                        % Extract temporal kernels
                        if params.regFlag == 0
                            counter = 2;
                        else
                            counter = 1;
                        end
                        for ne = 1:numel(params.eventNames)
                            if ~isnan(fits_real{na,ns}.betas.full)
                                clear thisKernel
                                thiscoeffs = fits_real{na,ns}.betas.full(counter:counter+numFun(ne)-1);
                                thisKernel = zscore(full(basisSets{ne}))*thiscoeffs;
                                eval(sprintf('temporalKernels.%s.%s(na,:,ns) = thisKernel'';',params.regions{nr},params.eventNames{ne}));
                            else
                                eval(sprintf('temporalKernels.%s.%s(na,:,ns) = nan.*zscore(full(basisSets{ne}))*ones(numel(counter:counter+numFun(ne)-1),1);',params.regions{nr},params.eventNames{ne}));
                            end
                            counter = counter+numFun(ne);
                        end
                    catch
                        if ns ~=4 && ~contains(thisIDs(na),'310-907')
                            keyboard
                        end
                        fprintf('\n\nNo data for mouse %s day %s \n\n',thisIDs{na},num2str(ns))
                        fits_real{na,ns}.betas.full       = NaN;
                        modelStats_real{na,ns}.corr.full = NaN;
                        modelStats_real{na,ns}.corr_noRefit.full = NaN;
                        modelStats_real{na,ns}.pvals = NaN;
                        modelStats_real{na,ns}.rsquared.full     = NaN;
                        modelStats_real{na,ns}.rsquared_adj.full = NaN;
                    end


                end
            end


            % Save data for each region
            saveName = fullfile(params.fbasename,'plots',params.regions{nr},sprintf('%s_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag)));
            save(saveName,'fits_real','fits_shuff','modelStats_real','modelStats_shuff','params','params');

            clear fits_real fits_shuff modelStats_real modelStats_shuff
        end

        % Save kernels and average fluorescence traces
        saveName = fullfile(params.fbasename,'plots',sprintf('%s_kernels_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag)));
        save(saveName,"meanFluor","temporalKernels","params");


        %% Model 2: Fit regression by session, concatenating across animals (option for mixed effects if regFlag == 0)
    case "model2"

        for nr = 1:numel(params.regions)
            thisIDs = unique(recs.ID(contains(recs.Region,params.regions{nr})));
            thisIDs(contains(thisIDs,excl.Var1)) = [];
            for ns = params.sessIDs
                lastwarn(''); % reset warning indicator
                %%% Initialize predictor matrix and response vector for original and shuffled data and animal and trial indicators for concatenating
                X = [];
                Y = [];
                aID = [];
                trialIndicator = [];

                counter = 0; % trial counter for cross validation
                for na = 1:numel(thisIDs)
                    %%% Concatenate predictor matrix, response variable, trial indicator, animal indicator across animals
                    dataLoc = fullfile(params.fbasename,params.regions{nr},thisIDs{na});
                    try % Catch statement for 310-907 who doesn't have data on session 4
                        %%% Create predictor matrix and response vector and concatenate
                        [thisX,thisY,thisTrialIndicator,numFun,basisSets,thisFluor,allFluor] = createPredictorMatrix_activeAvoid(params,params.eventNames,params.timeBack,params.timeForward,ns,dataLoc,0);
                        X                  = cat(1,X,thisX);
                        Y                  = cat(1,Y,thisY);
                        aID                = cat(1,aID,ones(size(thisY)).*na);
                        trialIndicator = cat(1,trialIndicator,thisTrialIndicator+counter);
                        counter        = counter+max(thisTrialIndicator);
                        %%% Extract avg fluorescence for each event
                        for ne = 1:numel(params.eventNames)
                            eval(sprintf('meanFluor.%s.%s(na,:,ns) = thisFluor{ne};',params.regions{nr},params.eventNames{ne}));
                        end
                    catch % 310-907 no day 4 data
                        if ns ~=4 && ~contains(thisIDs(na),'310-907')
                            keyboard
                        end
                    end
                end


                fprintf('Fitting %s day %s\n', params.regions{nr},  num2str(ns));

                %%% If using regularization, find lambda
                if params.regFlag > 0
                    fprintf('Finding lambda value \n')
                    lam=calcLambda(params.regFlag,X,Y,aID,trialIndicator);
                else
                    lam = NaN;
                end

                %%% Fit full and reduced models to real data.  If fitting without regularization, calculate F statistic for comparison between full and all reduced models (reduced models exclude all predictors associated w/ a particular event)
                [fits_real{1,ns},modelStats_real{1,ns},errFlag] = fitModel_activeAvoid(X,Y,params.eventNames,numFun,params.regFlag,aID,params.MEFlag,lam);

                if errFlag
                    % If fitting throws a warning, exclude fit (everything to NaN)
                    modelStats_real{1,ns}.corr.full = NaN;
                    modelStats_real{1,ns}.corr_noRefit.full = NaN;
                else
                    % Calculate correlation coefficient between real and estimated data for full and reduced models
                    if params.regFlag == 0
                        thisX = cat(2,ones(size(X,1),1),nanzscore(X,1)); % add column of ones to predictor matrix (intercept)
                        modelStats_real{1,ns} = calcModel_corr(thisX,nanzscore(Y,1),fits_real{1,ns},numFun, modelStats_real{1,ns},params.eventNames,params.regFlag,params.MEFlag,aID);
                    else
                        modelStats_real{1,ns} = calcModel_corr(nanzscore(X,1),nanzscore(Y,1),fits_real{1,ns},numFun, modelStats_real{1,ns},params.eventNames,params.regFlag,params.MEFlag,aID);
                    end
                end

                %%% Bootstrap fits to generate confidence intervals
                trials = unique(trialIndicator);
                if params.numShuff > 0
                    %%% Fit resampled data
                    fprintf('Bootstrapping');
                    if parFlag

                        parfor nss = 1:params.numShuff
                            %%% Generate resampled response vectors and predictor matrix
                            thisTrials = datasample(trials,numel(trials),1,"Replace",true); % resample trials with replacement
                            Y_shuff = cell2mat(arrayfun(@(x) Y(trialIndicator==x), thisTrials,'UniformOutput',false));
                            X_shuff = cell2mat(arrayfun(@(x) X(trialIndicator==x,:), thisTrials,'UniformOutput',false));
                            aID_shuff= cell2mat(arrayfun(@(x) aID(trialIndicator==x),thisTrials,'UniformOutput',false));
                            fprintf('.')
                            [tempFits_shuff{nss},tempStats_shuff{nss},shuffErr]  = fitModel_activeAvoid(X_shuff,Y_shuff,params.eventNames,numFun,params.regFlag,aID_shuff,params.MEFlag,lam);
                            if ~shuffErr
                                if params.regFlag == 0
                                    thisX = cat(2,ones(size(X_shuff,1),1),nanzscore(X_shuff,1));
                                    tempStats_shuff{nss} = calcModel_corr(thisX,nanzscore(Y_shuff,1),tempFits_shuff{nss},numFun, tempStats_shuff{nss},params.eventNames,params.regFlag,params.MEFlag,aID_shuff);
                                else
                                    tempStats_shuff{nss} = calcModel_corr(nanzscore(X_shuff,1),nanzscore(Y_shuff,1),tempFits_shuff{nss},numFun, tempStats_shuff{nss},params.eventNames,params.regFlag,params.MEFlag,aID_shuff);
                                end
                            else
                                tempStats_shuff{nss}.corr.full = NaN;
                            end
                        end
                        fits_shuff{1,ns} = tempFits_shuff;
                        modelStats_shuff{1,ns} = tempStats_shuff;
                        fprintf('\n')

                    else
                        for nss = 1:params.numShuff
                            %%% Generate resampled response vectors and predictor matrix
                            thisTrials = datasample(trials,numel(trials),1,"Replace",true); % resample trials with replacement
                            Y_shuff  = cell2mat(arrayfun(@(x) Y(trialIndicator==x), thisTrials,'UniformOutput',false));
                            X_shuff  = cell2mat(arrayfun(@(x) X(trialIndicator==x,:), thisTrials,'UniformOutput',false));
                            aID_shuff= cell2mat(arrayfun(@(x) aID(trialIndicator==x),thisTrials,'UniformOutput',false));
                            fprintf('.')
                            [fits_shuff{1,ns}{nss},modelStats_shuff{1,ns}{nss},shuffErr]  = fitModel_activeAvoid(X_shuff,Y_shuff,params.eventNames,numFun,params.regFlag,aID_shuff,params.MEFlag,lam);
                            if ~shuffErr
                                if params.regFlag == 0
                                    thisX = cat(2,ones(size(X_shuff,1),1),nanzscore(X_shuff,1));
                                    modelStats_shuff{1,ns}{nss} = calcModel_corr(thisX,nanzscore(Y_shuff,1),fits_shuff{1,ns}{nss},numFun, modelStats_shuff{1,ns}{nss},params.eventNames,params.regFlag,params.MEFlag,aID_shuff);
                                else
                                    thisX = nanzscore(X_shuff,1);
                                    modelStats_shuff{1,ns}{nss} = calcModel_corr(thisX,nanzscore(Y_shuff,1),fits_shuff{1,ns}{nss},numFun, modelStats_shuff{1,ns}{nss},params.eventNames,params.regFlag,params.MEFlag,aID_shuff);
                                end
                            else
                                modelStats_shuff{1,ns}{nss}.corr.full = NaN;
                            end
                        end
                        fprintf('\n')
                    end
                else
                    fits_shuff{1,ns}{1} = NaN;
                    modelStats_shuff{1,ns}{1} = NaN;
                end



                %%% Organize data for plotting
                % Extract temporal kernels
                if params.regFlag == 0
                    counter = 2;
                else
                    counter = 1;
                end
                for ne = 1:numel(params.eventNames)
                    if ~isnan(fits_real{1,ns}.betas.full)
                        clear thisKernel
                        thiscoeffs = fits_real{1,ns}.betas.full(counter:counter+numFun(ne)-1);
                        thisKernel = zscore(full(basisSets{ne}))*thiscoeffs;
                        eval(sprintf('temporalKernels.%s.%s(1,:,ns) = thisKernel'';',params.regions{nr},params.eventNames{ne}));
                        if params.numShuff > 0
                            for nss = 1:params.numShuff
                                thiscoeffs = fits_shuff{1,ns}{nss}.betas.full(counter:counter+numFun(ne)-1);
                                thisKernel = zscore(full(basisSets{ne}))*thiscoeffs;
                                eval(sprintf('shuffKernels.%s.%s(nss,:,ns) = thisKernel'';',params.regions{nr},params.eventNames{ne}));
                            end
                        end
                    else
                        eval(sprintf('temporalKernels.%s.%s(1,:,ns) = nan.*zscore(full(basisSets{ne}))*ones(numel(counter:counter+numFun(ne)-1),1);',params.regions{nr},params.eventNames{ne}));
                    end
                    counter = counter+numFun(ne);
                end
            end
        
        % Save data for each region
        saveName = fullfile(params.fbasename,'plots',params.regions{nr},sprintf('%s_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag)));
        save(saveName,'fits_real','fits_shuff','modelStats_real','modelStats_shuff','params','params');

        clear fits_real fits_shuff modelStats_real modelStats_shuff
        end
        % Save kernels and average fluorescence traces
        saveName = fullfile(params.fbasename,'plots',sprintf('%s_kernels_numShuff%s_Fs%s_ME%s_reg%s.mat',params.model,num2str(params.numShuff),num2str(params.newFs),num2str(params.MEFlag),num2str(params.regFlag)));
        save(saveName,"meanFluor","temporalKernels","shuffKernels","params");
end
end



%% Functions

%% Find best lambda using full data set based on cross-validated correlation between real and estimated data
function lam = calcLambda(regFlag,X,Y,aID,trialIndicator)

% Set up partitions for determining lambda
% Find trials for each animal
IDs = unique(aID);
counter = 1;
trials = [];
group = [];
for na = 1:numel(IDs)
    trainTrials = unique(trialIndicator(aID==IDs(na)));
    trials = cat(1,trials,trainTrials+counter-1);
    trialIndicator(aID==IDs(na)) = trialIndicator(aID==IDs(na))+counter-1;
    group = cat(1,group,ones(size(trainTrials)).*na);
    counter = counter+max(trainTrials);
end

cv = cvpartition(group,"KFold",5,"Stratify",true); % Create cross-validation partitions with stratification by animal or day (if concatenated)

if regFlag == 1
    lambdas = 0:100:10000; % what lambda values to test for ridge regression
    % find lambda that gives best cross-validated fit
    errVec = zeros(1,cv.NumTestSets);
    for nf = 1:cv.NumTestSets
        thisTrain = find(training(cv,nf));
        thisTest  = find(test(cv,nf));
        thisY = zscore(Y(sum(trialIndicator==thisTrain',2)==1),[],1);
        thisX = zscore(X(sum(trialIndicator==thisTrain',2)==1,:),[],1);
        B = ridge(thisY,thisX,lambdas);
        for nb = 1:size(B,2)
            yhat = zscore(X(sum(trialIndicator==thisTest',2)==1,:),[],1)*B(:,nb); % estimated data
            testY = zscore(Y(sum(trialIndicator==thisTest',2)==1,:),[],1);
            lamMSE(nb,nf) = mean(sum((testY-yhat).^2));
        end
        [~,lam(nf)] = min(lamMSE(:,nf));
        lam(nf) = lambdas(lam(nf));
        if ~isempty(lastwarn)
            errVec(nf) = 1;
            lastwarn('');
        end
    end
    lam = mean(lam(~errVec));


elseif regFlag == 2
    % Find lambda that gives best cross-validated fit
    %[B,etc] = lasso(zscore(X,[],1),zscore(Y,[],1),'NumLambda',1000,'CV',10);
    %lam = etc.Lambda(etc.Index1SE);
    errVec = zeros(1,cv.NumTestSets);
    for nf = 1:cv.NumTestSets
        thisTrain = find(training(cv,nf));
        thisTest  = find(test(cv,nf));
        thisY = zscore(Y(sum(trialIndicator==thisTrain',2)==1),[],1);
        thisX = zscore(X(sum(trialIndicator==thisTrain',2)==1,:),[],1);
        [B,etc] = lasso(thisX,thisY,'NumLambda',1000);
        for nb = 1:size(B,2)
            yhat = zscore(X(sum(trialIndicator==thisTest',2)==1,:),[],1)*B(:,nb); % estimated data
            testY = zscore(Y(sum(trialIndicator==thisTest',2)==1,:),[],1);
            lamMSE(nb,nf) = mean(sum((testY-yhat).^2));
            allLam(nb,nf) = etc.Lambda(nb);
        end
        [~,lam(nf)] = min(lamMSE(:,nf));
        lam(nf) = etc.Lambda(lam(nf));
        if ~isempty(lastwarn)
            errVec(nf) = 1;
            lastwarn('');
        end
    end
    lam = mean(lam(~errVec));
end

end




