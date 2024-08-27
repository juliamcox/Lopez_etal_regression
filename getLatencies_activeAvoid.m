function params = getLatencies_activeAvoid(params,fbasename)

%% Find animal IDs
% Animal IDs
% Load list of recordings
recs = readtable(fullfile(fbasename,'Photometry_Identification.xlsx'));
recs.Properties.VariableNames = {'ID';'Fname';'Day';'Region'};
Fnames = recs.Fname;
Fnames = cellfun(@(x) x(2:end-1), Fnames, 'UniformOutput', false);
recs.Fname = string(Fnames);
% Excluded mice
excl = readtable(fullfile(fbasename,'ExcludedMice_MouseID.csv'),ReadVariableNames=false);


for nr = 1:numel(params.regions)
    thisIDs = unique(recs.ID(contains(recs.Region,params.regions{nr})));
    thisIDs(contains(thisIDs,excl.Var1)) = [];
    avoidLatency = cell(1,7);
    escapeLatency= cell(1,7);
    for na = 1:numel(thisIDs)
        load(fullfile(fbasename,params.regions{nr},thisIDs{na},'behavEvents.mat'),'behavEvents');
        for ns = 1:7
            thisEvent = [];
            thisType  = [];
            eventNames = behavEvents.EventName(behavEvents.Day==ns);
            for ne = 1:numel(eventNames)

                thisEvent = cat(1,thisEvent,behavEvents.TS{behavEvents.Day==ns&behavEvents.EventName==eventNames(ne)});
                thisType  = cat(1, thisType,repmat(eventNames(ne),size(behavEvents.TS{behavEvents.Day==ns&behavEvents.EventName==eventNames(ne)})));
            end
            % Sort events
            [thisEvent,idx] = sort(thisEvent);
            thisType        = thisType(idx);
            if ~isempty(idx)
                % Figure out which events belong to which trials
                temp = diff(thisEvent);
                tempIdx = cat(1,1,find(temp>=30)+1,numel(thisEvent)+1);
                trialIndicator = cell2mat(arrayfun(@(x,y,z) ones(size(x:y)).*z,tempIdx(1:end-1),tempIdx(2:end)-1,[1:numel(tempIdx)-1]','UniformOutput',false)')';

                thisAvoidTrials = trialIndicator(thisType=="CueAvoid");
                try
                    thisAvoid       = arrayfun(@(x) thisEvent(trialIndicator==x&thisType=="AvoidCross") - thisEvent(trialIndicator==x&thisType=="CueAvoid"), thisAvoidTrials);
                catch
                    if numel(unique(thisAvoidTrials)) ~= numel(thisAvoidTrials)
                        % ignore it for now
                        [~,i] = unique(thisAvoidTrials);
                        dupTrial = thisAvoidTrials(find(find(diff(i)>1)));
                        thisAvoid       = arrayfun(@(x) thisEvent(find(trialIndicator==x&thisType=="AvoidCross",1,'first')) - thisEvent(find(trialIndicator==x&thisType=="CueAvoid",1,'first')), thisAvoidTrials);
                        thisAvoid(dupTrial) = [];
                    else
                        keyboard
                    end
                end
                avoidLatency{ns}=  cat(1,avoidLatency{ns}, thisAvoid);

                thisEscapeTrials = trialIndicator(thisType=="Shock");
                try
                    thisEscape       = arrayfun(@(x) thisEvent(trialIndicator==x&thisType=="EscapeCross") - thisEvent(trialIndicator==x&thisType=="Shock"), thisEscapeTrials);
                catch
                    thisEscape       = cell2mat(arrayfun(@(x) thisEvent(trialIndicator==x&thisType=="EscapeCross") - thisEvent(trialIndicator==x&thisType=="Shock"), thisEscapeTrials,'UniformOutput',false));

                    
                end
                escapeLatency{ns}=  cat(1,escapeLatency{ns}, thisEscape);

            end
        end

    end
    eval(sprintf('params.escapeLatency.%s = escapeLatency;',params.regions{nr}))
    eval(sprintf('params.avoidLatency.%s = avoidLatency;',params.regions{nr}));
end


        
        
        
        
        
        
        
        
        
        
        
        
        
        