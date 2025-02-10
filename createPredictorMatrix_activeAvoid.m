function [X,Y,trialIndicator,numFun,basisSets,meanFluor,fluor] = createPredictorMatrix_activeAvoid(params,eventNames,timeBack,timeForward,sessIDs,dataLoc,shuffFlag,speedFlag)

% Creates predictor matrix and response vector (not z-scored) 
% predictor matrix is made using a b-spline basis set 

%%% Inputs %%%
% params: structure containing original sampling rate (params.Fs), sampling rate to downsample to (params.newFs)
% eventNames: cell array of event names to include in the regression 
% timeBack: time before event (in sec)
% timeForward: time after the event (in sec)
% sessIDs: which days to include
% dataLoc: location of data 
% shuffFlag: circularly shift the response vector? 0 or 1 
% speedFlag: is speed a predictor? 

%%% Outputs %%%
% X: predictor matrix
% Y: response vector 
% trialIndicator: indicates which trial each data point belongs to (for trial-wise coross validation). if numel(sessIDs)>1 trial count is continuous
% numFun: number of functions used for each predictor (params.numBasis per second) 
% basisSets: the basis sets used to create the predictor matrix
% meanFluor: average event triggered fluorescence for each event in eventNames from timeBack:timeForward
% fluor: trial by trial event triggered fluorescence for each event in eventNames from timeBack:timeForward

if nargin < 8
    speedFlag = false;
end


%% Extract response vector and predictor matrix, concatenating across sessions in sessIDs 
if speedFlag
    % Initialize predictor matrix and response vector
    X = [];
    Y = [];
    trialIndicator = []; % which trial is which for cross validation

    fluor = cell(size(eventNames));
    meanFluor = cell(size(eventNames));

    counter = 0;
    for ns = 1:numel(sessIDs)

        thisX = []; % initialize predictor matrix
        idx   = []; % initialize list of in-trial indices

        %%% Load photometry trace
        load(fullfile(dataLoc,sprintf('DAFluor_day%s.mat',num2str(sessIDs(ns)))));
        %%% Downsample photometry trace if needed
        if params.Fs ~= params.newFs
            thisY    = resample(dataTable.dff,params.newFs,params.Fs);
            thisSpeed = resample(dataTable.Speed,params.newFs,params.Fs);
        else
            thisY     = dataTable.dff;
            thisSpeed = dataTable.Speed;
        end
        thisAccel = cat(1,NaN,diff(thisSpeed)./(1/params.newFs));
        %%% Circularly shift data if necessary
        if shuffFlag
            thisY = circshift(thisY,randi(size(thisY,1)-params.Fs));
        end
        for ne = 1:numel(eventNames)
            try
                thisEvents = dataTable.Time(eval(sprintf('dataTable.%s',eventNames{ne}))==1);
                thisEvents = round(thisEvents.*params.newFs); % find frame number for each event

                %%% Extract mean fluorescence
                clear thisFluor
                for nt = 1:numel(thisEvents)
                    if thisEvents(nt)+(timeForward(ne)*params.newFs) > size(thisY,1)
                        thisFluor(nt,:) = cat(2,thisY(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):end)', nan(1,numel(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):round(thisEvents(nt)+(timeForward(ne)*params.newFs)))- numel(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):size(thisY,1))));
                        idx = cat(2,idx,(thisEvents(nt)-timeBack(ne)*params.newFs):size(thisY,1));
                    else
                        thisFluor(nt,:) = thisY(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):round(thisEvents(nt)+(timeForward(ne)*params.newFs)));
                        idx = cat(2,idx,thisEvents(nt)-(timeBack(ne)*params.newFs):thisEvents(nt)+(timeForward(ne)*params.newFs));
                    end
                end
                fluor{ne} = cat(1,fluor{ne},thisFluor);
                meanFluor{ne} = cat(1,meanFluor{ne},mean(thisFluor,1,'omitnan'));
            catch % if there were no trials for this event
                thisEvents = [];
                if sum(behavEvents.Day==sessIDs(ns)&behavEvents.EventName==eventNames{ne}) ~= 0
                    keyboard
                end

                %%% Extract mean fluorescence
                clear thisFluor
                thisFluor = nan(1,numel(-(timeBack(ne)*params.newFs):+(timeForward(ne)*params.newFs)));
                fluor{ne} = cat(1,fluor{ne},thisFluor);
                meanFluor{ne} = cat(1,meanFluor{ne},mean(thisFluor,1,'omitnan'));

            end

            %%% Generate basis set
            range  = [-timeBack(ne)*params.newFs timeForward(ne)*params.newFs];
            numFun(ne) = round(params.numBasis*numel(range(1):range(2))./params.newFs);
            order  = 4;
            basisobj = create_bspline_basis(range,numFun(ne),order);
            bs = getbasismatrix(range(1):range(2),basisobj);
            basisSets{ne} = bs;

            %%% Create predictor matrix
            % Create binary vector of event start times
            temp = zeros(size(thisY));
            temp(thisEvents-round(timeBack(ne)*params.newFs)) = 1;

            % Convolve with the basis set
            x=[];
            for nb = 1:numFun(ne)
                thisPred = conv(temp,full(bs(:,nb)));
                x = cat(2,x,thisPred(1:size(temp,1)));
            end
            %%% Concatenate
            thisX = cat(2,thisX,x);

        end
        if params.model == "model5"
            idx = 1:size(thisX,1);
        else
            idx = unique(idx);
            idx = round(idx);
        end
        thisY = thisY(idx);
        thisX = thisX(idx,:);
        thisX = cat(2,thisX,thisSpeed(idx));
        if params.model == "model4" || params.model == "model5"
            thisX = cat(2,thisX,thisAccel(idx));
        end


        %%% Create trial inidicator
        temp = diff(idx);
        tempIdx = cat(2,1,find(temp>=(30 - max(timeBack(contains(eventNames,'Cue'))) - max(timeForward(contains(eventNames,'Cross'))))*params.newFs)+1,numel(idx)+1);

        thisIndicator = cell2mat(arrayfun(@(x,y,z) ones(size(x:y))'.*z,tempIdx(1:end-1),tempIdx(2:end)-1,1:numel(tempIdx)-1,'UniformOutput',false)');
        trialIndicator= cat(1,trialIndicator,thisIndicator+counter);
        counter = counter+max(trialIndicator);


        %%% Concatenate X and Y
        Y = cat(1,Y,thisY);
        X = cat(1,X,thisX);

    end

    for ne = 1:numel(eventNames)
        meanFluor{ne} = mean(meanFluor{ne},1,'omitnan');
    end


else
    % Initialize predictor matrix and response vector
    X = [];
    Y = [];
    trialIndicator = []; % which trial is which for cross validation

    fluor = cell(size(eventNames));
    meanFluor = cell(size(eventNames));

    % Load timestamps
    load(fullfile(dataLoc, 'behavEvents.mat'),'behavEvents');

    counter = 0;
    for ns = 1:numel(sessIDs)

        thisX = []; % initialize predictor matrix
        idx   = []; % initialize list of in-trial indices

        %%% Load photometry trace
        load(fullfile(dataLoc,sprintf('DAFluor_day%s.mat',num2str(sessIDs(ns)))));
        %%% Downsample photometry trace if needed
        if params.Fs ~= params.newFs
            thisY   = resample(data,params.newFs,params.Fs);
        else
            thisY    = data;
        end
        %%% Circularly shift data if necessary
        if shuffFlag
            thisY = circshift(thisY,randi(size(thisY,1)-params.Fs));
        end

        for ne = 1:numel(eventNames)
            try
                thisEvents = behavEvents.TS{behavEvents.Day==sessIDs(ns)&behavEvents.EventName==eventNames{ne}};
                thisEvents = round(thisEvents.*params.newFs); % find frame number for each event

                %%% Extract mean fluorescence
                clear thisFluor
                for nt = 1:numel(thisEvents)
                    if thisEvents(nt)+(timeForward(ne)*params.newFs) > size(thisY,1)
                        thisFluor(nt,:) = cat(2,thisY(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):end)', nan(1,numel(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):round(thisEvents(nt)+(timeForward(ne)*params.newFs)))- numel(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):size(thisY,1))));
                        idx = cat(2,idx,(thisEvents(nt)-timeBack(ne)*params.newFs):size(thisY,1));
                    else
                        thisFluor(nt,:) = thisY(round(thisEvents(nt)-(timeBack(ne)*params.newFs)):round(thisEvents(nt)+(timeForward(ne)*params.newFs)));
                        idx = cat(2,idx,thisEvents(nt)-(timeBack(ne)*params.newFs):thisEvents(nt)+(timeForward(ne)*params.newFs));
                    end
                end
                fluor{ne} = cat(1,fluor{ne},thisFluor);
                meanFluor{ne} = cat(1,meanFluor{ne},mean(thisFluor,1,'omitnan'));
            catch % if there were no trials for this event
                thisEvents = [];
                if sum(behavEvents.Day==sessIDs(ns)&behavEvents.EventName==eventNames{ne}) ~= 0
                    keyboard
                end

                %%% Extract mean fluorescence
                clear thisFluor
                thisFluor = nan(1,numel(-(timeBack(ne)*params.newFs):+(timeForward(ne)*params.newFs)));
                fluor{ne} = cat(1,fluor{ne},thisFluor);
                meanFluor{ne} = cat(1,meanFluor{ne},mean(thisFluor,1,'omitnan'));

            end

            %%% Generate basis set
            range  = [-timeBack(ne)*params.newFs timeForward(ne)*params.newFs];
            numFun(ne) = round(params.numBasis*numel(range(1):range(2))./params.newFs);
            order  = 4;
            basisobj = create_bspline_basis(range,numFun(ne),order);
            bs = getbasismatrix(range(1):range(2),basisobj);
            basisSets{ne} = bs;

            %%% Create predictor matrix
            % Create binary vector of event start times
            temp = zeros(size(thisY));
            temp(thisEvents-round(timeBack(ne)*params.newFs)) = 1;

            % Convolve with the basis set
            x=[];
            for nb = 1:numFun(ne)
                thisPred = conv(temp,full(bs(:,nb)));
                x = cat(2,x,thisPred(1:size(temp,1)));
            end
            %%% Concatenate
            thisX = cat(2,thisX,x);

        end

        idx = unique(idx);
        idx = round(idx);
        thisY = thisY(idx);
        thisX = thisX(idx,:);


        %%% Create trial inidicator
        temp = diff(idx);
        tempIdx = cat(2,1,find(temp>=(30 - max(timeBack(contains(eventNames,'Cue'))) - max(timeForward(contains(eventNames,'Cross'))))*params.newFs)+1,numel(idx)+1);

        thisIndicator = cell2mat(arrayfun(@(x,y,z) ones(size(x:y))'.*z,tempIdx(1:end-1),tempIdx(2:end)-1,1:numel(tempIdx)-1,'UniformOutput',false)');
        trialIndicator= cat(1,trialIndicator,thisIndicator+counter);
        counter = counter+max(trialIndicator);


        %%% Concatenate X and Y
        Y = cat(1,Y,thisY);
        X = cat(1,X,thisX);

    end

    for ne = 1:numel(eventNames)
        meanFluor{ne} = mean(meanFluor{ne},1,'omitnan');
    end
end