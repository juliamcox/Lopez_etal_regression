function [fits,modelStats,errFlag] = fitModel_activeAvoid(X,Y,eventNames,numFun,regFlag,aIDs,MEFlag,lam,speedFlag)

% Fit (mixed-effects) regression with or without regularization

%%% Inputs:
% X: predictor matrix  
% Y: response vector
% numFun: number of basis set functions corresponding to each event 
% regFlag: fit with regularization? 0 = no, 1 = ridge, 2 = lasso 
% MEFlag: fit mixed effects model? Do not use regularization in this case 
% aIDs: vector the same size as Y indicating which animal each data point comes from 
% lam: regularization parameter 
% speedFlag: is there a speed predictor?

%%% Outputs:
% fits: structure containing coefficients for full and reduced models 
% modelStats: structure containing goodness of fit metrics (when regFlag == 0) for full and reduced models w/ and w/o refitting 
% errFlag: did the fit throw and error? 


errFlag = 0;
lastwarn(''); % reset warning indicator

%% Fit full and reduced models (without each event)  

if MEFlag
    %% fit mixed effects model with random intercept for animal ID
    % configure inputs
    M = array2table(nanzscore(X,1));
    vNames = M.Properties.VariableNames;
    M.aID = aIDs;
    M.Y = nanzscore(Y);
    f = 'Y ~ ';
    for nv = 1:numel(vNames)
        f = cat(2,f,sprintf('%s +',vNames{nv}));
    end
    % fit model with random intercept for subject
    f = cat(2,f,'(1|aID)');
    mdl = fitlme(M,f);
    modelStats.rsquared.full = mdl.Rsquared.Ordinary;
    modelStats.rsquared_adj.full = mdl.Rsquared.Adjusted;
    % Extract coefficients
    if isempty(lastwarn)
        fits.betas.full               = mdl.Coefficients.Estimate;
        fits.randomEffects.full = randomEffects(mdl);
    else
        % if the model fit threw a warning, exclude 
        fits.betas.full              = nan(size(X,2),1);
        modelStats.rsquared.full     = NaN;
        modelStats.rsquared_adj.full = NaN;
        lastwarn(''); % reset warning indicator
        % If there is a warning (i.e., matrix not full rank) don't run the reduced models/shuffles
        errFlag = 1;
        return
    end
    % fit reduced models 
    counter = 1;
    for ne = 1:numel(eventNames)
        % Exclude predictors for this event
        thisPreds = counter:counter+numFun(ne)-1;
        % remove predictors for event
        thisX = X;
        thisX(:,thisPreds) = [];
        eval(sprintf('modelStats.Fstat.%s = NaN;',eventNames{ne})); % no f-statistic for mixed effects model 
        M = array2table(nanzscore(thisX,1));
        vNames = M.Properties.VariableNames;
        M.aID = aIDs;
        M.Y = nanzscore(Y);
        f = 'Y ~ ';
        for nv = 1:numel(vNames)
            f = cat(2,f,sprintf('%s +',vNames{nv}));
        end
        f = cat(2,f,'(1|aID)');
        mdl = fitlme(M,f);
        eval(sprintf('fits.betas.%s = mdl.Coefficients.Estimate;',eventNames{ne}));
        eval(sprintf('fits.randomEffects.%s = randomEffects(mdl);',eventNames{ne}));
        counter = counter+numFun(ne);
        eval(sprintf('modelStats.rsquared.%s = mdl.Rsquared.Ordinary;',eventNames{ne}));
        eval(sprintf('modelStats.rsquared_adj.%s = mdl.Rsquared.Adjusted;',eventNames{ne}));
    end



elseif regFlag == 0
    %% Fit linear regression without regularization
    mdl   = fitlm(nanzscore(X,1),nanzscore(Y,1));
    modelStats.rsquared.full = mdl.Rsquared.Ordinary;
    modelStats.rsquared_adj.full = mdl.Rsquared.Adjusted;
    % Extract coefficients
    if isempty(lastwarn)
        fits.betas.full = mdl.Coefficients.Estimate;
    else
        % if the model fit threw a warning, exclude
        fits.betas.full              = nan(size(X,2)+1,1);
        modelStats.rsquared.full     = NaN;
        modelStats.rsquared_adj.full = NaN;
        lastwarn(''); % reset warning indicator
        % If there is a warning (i.e., matrix not full rank) don't run the reduced models/shuffles
        errFlag = 1;
    end
    counter = 1;
    for ne = 1:numel(eventNames)
        % Exclude predictors for this event
        thisPreds = counter:counter+numFun(ne)-1;
        % remove predictors for event
        thisX = X;
        thisX(:,thisPreds) = [];
        % Fit reduced model
        if errFlag
            eval(sprintf('fits.betas.%s = nan(size(thisX,2)+1,1);',eventNames{ne}));
        else
            mdl   = fitlm(nanzscore(thisX,1),nanzscore(Y,1));
            eval(sprintf('modelStats.rsquared.%s = mdl.Rsquared.Ordinary;',eventNames{ne}));
            eval(sprintf('modelStats.rsquared_adj.%s = mdl.Rsquared.Adjusted;',eventNames{ne}));
            eval(sprintf('fits.betas.%s = mdl.Coefficients.Estimate;',eventNames{ne}));

            % Calculate F statistic comparing full and reduced model
            [Fstat, Fstat_noRefit] = get_fStat(nanzscore(X,1),nanzscore(thisX,1),nanzscore(Y,1),fits.betas.full,mdl.Coefficients.Estimate,thisPreds);
            eval(sprintf('modelStats.Fstat.%s = Fstat;',eventNames{ne}));
            eval(sprintf('modelStats.Fstat_noRefit.%s = Fstat_noRefit;',eventNames{ne}));
        end
        counter = counter+numFun(ne);
    end


%% Fit linear regression with regularization
else
    modelStats = [];
    modelStats.lam = lam;
%% ridge regularization
    if regFlag == 1
        
        % Fit the full model
       fits.betas.full = ridge(nanzscore(Y,1),nanzscore(X,1),lam);
       if ~isempty(lastwarn)
           keyboard
       end

        % Fit the reduced models
        counter = 1;
        for ne = 1:numel(eventNames)
            % Exclude predictors for this event
            thisPreds = counter:counter+numFun(ne)-1;
            % remove predictors for event
            thisX = X;
            thisX(:,thisPreds) = [];
            % Fit reduced model
            B   = ridge(nanzscore(Y,1),nanzscore(thisX,1),lam);
            eval(sprintf('fits.betas.%s = B;',eventNames{ne}));
            counter = counter+numFun(ne);
        end

%% Lasso regularization 
    elseif regFlag == 2 
         % fit the full model
         try
         fits.betas.full = lasso(nanzscore(X,1),nanzscore(Y,1),'Lambda',lam);
         catch
             keyboard
         end
         if ~isempty(lastwarn)
             keyboard
         end
         % Fit the reduced models
         counter = 1;
         for ne = 1:numel(eventNames)
             % Exclude predictors for this event
             thisPreds = counter:counter+numFun(ne)-1;
             % remove predictors for event
             thisX = X;
             thisX(:,thisPreds) = [];
             % Fit reduced model
             B   = lasso(nanzscore(thisX,1),nanzscore(Y,1),'Lambda',lam);
             eval(sprintf('fits.betas.%s = B;',eventNames{ne}));
             counter = counter+numFun(ne);
         end

         if speedFlag
             for ne = counter:size(X,2)
                thisX = X;
                thisX(:,ne) = [];
                B   = lasso(nanzscore(thisX,1),nanzscore(Y,1),'Lambda',lam);
                eval(sprintf('fits.betas.speed%s = B;',num2str(counter-size(X,2))))
             end

         end

    end
end

end



