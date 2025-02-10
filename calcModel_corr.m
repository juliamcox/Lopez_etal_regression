function modelStats = calcModel_corr(X,Y,fits,numFun, modelStats, eventNames,regFlag,MEFlag,aIDs,speedFlag)



%% Calculate correlation coefficients between real and estimated data 

% calculate corr coeff for full model
if MEFlag
    g = unique(aIDs)';
    Z = zeros(size(Y,1),numel(g));
    for na = 1:numel(g)
        Z(aIDs==na,na) = 1;
    end
    yhat = X*fits.betas.full + Z*fits.randomEffects.full;
else
    yhat = X*fits.betas.full;
end
modelStats.corr.full = corr(yhat,Y,'rows','complete');


% Calculate correlation coefficient for the reduced models
if regFlag == 0
    counter = 2;
else
    % if fit with regularization, no intercept 
    counter = 1;
end
for ne = 1:numel(eventNames)
    % Exclude predictors for this event
    thisX = X;
    thisX(:,counter:counter+numFun(ne)-1) = [];
    thisB = eval(sprintf('fits.betas.%s;',eventNames{ne}));
    if MEFlag
        thisRE = eval(sprintf('fits.randomEffects.%s;',eventNames{ne}));
        yhat = thisX*thisB + Z*thisRE;
    else
        yhat  = thisX*thisB; % estimated data
    end
    eval(sprintf('modelStats.corr.%s = corr(yhat,Y,''rows'',''complete'');',eventNames{ne}));
    counter = counter+numFun(ne);
end

if speedFlag
    for ne = counter:size(X,2)
        % Exclude predictors for this event
        thisX = X;
        thisX(:,counter) = [];
        thisB = eval(sprintf('fits.betas.speed%s;',num2str(size(X,2)-counter)));
        if MEFlag
            thisRE = eval(sprintf('fits.randomEffects.speed%s;',num2str(size(X,2)-counter)));
            yhat = thisX*thisB + Z*thisRE;
        else
            yhat  = thisX*thisB; % estimated data
        end
        eval(sprintf('modelStats.corr.speed%s = corr(yhat,Y,''rows'',''complete'');',num2str(size(X,2)-counter)));
    end
end

% Calculate correlation coefficient for the reduced models without refitting 
if regFlag == 0
    counter = 2;
else
    % if fit with regularization, no intercept 
    counter = 1;
end
for ne = 1:numel(eventNames)
    % Exclude predictors for this event
    thisX = X;
    thisX(:,counter:counter+numFun(ne)-1) = [];
    thisB = fits.betas.full;
    thisB(counter:counter+numFun(ne)-1) = []; 
    if MEFlag
        thisRE = fits.randomEffects.full;
        yhat = thisX*thisB + Z*thisRE;
    else
        yhat  = thisX*thisB; % estimated data
    end
    eval(sprintf('modelStats.corr_noRefit.%s = corr(yhat,Y,''rows'',''complete'');',eventNames{ne}));
    counter = counter+numFun(ne);
end

if speedFlag
    for ne = counter:size(X,2)
        % Exclude predictors for this event
        thisX = X;
        thisX(:,counter) = [];
        thisB = fits.betas.full;
        thisB(counter) = [];
        if MEFlag
            thisRE = fits.randomEffects.full;
            yhat = thisX*thisB + Z*thisRE;
        else
            yhat  = thisX*thisB; % estimated data
        end
        eval(sprintf('modelStats.corr_noRefit.speed%s = corr(yhat,Y,''rows'',''complete'');',num2str(size(X,2)-counter)));
        counter = counter+1;
    end
end

