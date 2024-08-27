function modelStats = calcModel_corr(X,Y,fits,numFun, modelStats, eventNames,regFlag,MEFlag,aIDs)



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
modelStats.corr.full = corr(yhat,Y);


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
    eval(sprintf('modelStats.corr.%s = corr(yhat,Y);',eventNames{ne}));
    counter = counter+numFun(ne);
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
    eval(sprintf('modelStats.corr_noRefit.%s = corr(yhat,Y);',eventNames{ne}));
    counter = counter+numFun(ne);
end

