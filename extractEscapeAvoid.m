function extractEscapeAvoid(fbasename,saveLoc,params)

% Location of photometry traces
FPLoc   = fullfile(fbasename,'FP_forJulia');
% Excluded mice 
excl = readtable(fullfile(fbasename,'ExcludedMice_MouseID.csv'),ReadVariableNames=false);

% Make full list of file locations
allFiles = [];
for nr = 1:numel(params.regions)
    flist = dir(fullfile(FPLoc,params.regions{nr}));
    flist = {flist(:).name};
    flist = flist(contains(flist,'Cohort'));
    for nf = 1:numel(flist)
        thisF = dir(fullfile(FPLoc,params.regions{nr},flist{nf}));
        thisF = thisF(contains({thisF(:).name},'Day'));
        for nff = 1:numel(thisF)
            thisFF = dir(fullfile(FPLoc,params.regions{nr},flist{nf},thisF(nff).name));
            thisFF = {thisFF(:).name};
            temp = cellfun(@(x) fullfile(FPLoc,params.regions{nr},flist{nf},thisF(nff).name,x),thisFF,'UniformOutput',false)';
            allFiles = cat(1,allFiles,temp);
        end
    end
end

allFiles = string(allFiles);



% Load list of recordings
recs = readtable(fullfile(FPLoc,'Photometry_Identification.xlsx'));
recs.Properties.VariableNames = {'ID';'Fname';'Day';'Region'};
Fnames = recs.Fname;
for nf = 1:numel(Fnames)
    idx = strfind(Fnames{nf},'''');
    Fnames{nf}(idx) = [];
end
recs.Fname = string(Fnames);



for nr = 1:numel(params.regions)
    IDs  = unique(recs.ID(recs.Region==string(params.regions{nr})));
    IDs(contains(IDs,excl.Var1)) = [];
    for na = 1:numel(IDs)
        behavEvents = array2table(zeros(0,3), 'VariableNames',{'TS';'Day';'EventName'}); % initialize behavior events table 
        fnames = recs.Fname(contains(recs.ID,IDs{na})); % find file locations for each animal
        days   = recs.Day(contains(recs.ID,IDs{na})); % find recording day for each animal
        % sort by session
        [~,idx]= sort(days);
        fnames = fnames(idx);
        days   = days(idx);

        for nf = 1:numel(fnames)
            % Extract photometry trace
            thisName = allFiles(contains(allFiles,fnames{nf}));
            if numel(thisName) > 1
                if contains(lower(fnames{nf}),'copy')
                    keyboard
                    thisName = thisName(contains(lower(thisName),'copy'));
                else
                    thisName = thisName(~contains(lower(thisName),'copy'));
                end
            end
            try
            floc = dir(thisName);
            catch
                keyboard
            end
            floc = string({floc(:).name});
            floc = floc(contains(floc,'output_1'));
            data = h5read(fullfile(thisName,floc,'dff_DA.hdf5'),'/data');
            % save mat file
            if ~isdir(fullfile(saveLoc,params.regions{nr},IDs{na}))
                mkdir(fullfile(saveLoc,params.regions{nr},IDs{na}))
            end
            save(fullfile(saveLoc,params.regions{nr},IDs{na},sprintf('DAFluor_day%s.mat',num2str(days(nf)))),'data');

            % Extract timestamps


            for ne = 1:numel(params.eventNames_all)
                flist = dir(fullfile(fbasename,params.regions{nr},params.eventNames_all{ne},sprintf('Day%s',num2str(days(nf)))));
                flist = string({flist(:).name});
                flist = flist(contains(flist,IDs{na}));
                try
                    variableNames = who('-file',fullfile(fbasename,params.regions{nr},params.eventNames_all{ne},sprintf('Day%s',num2str(days(nf))),flist));

                    temp = table(string(params.eventNames_all{ne}),'VariableNames',{'EventName'});
                    temp.Day = days(nf);
                    try
                        load(fullfile(fbasename,params.regions{nr},params.eventNames_all{ne},sprintf('Day%s',num2str(days(nf))),flist),variableNames{contains(variableNames,'ts_matrix')});
                        temp.TS = {ts_matrix};
                        clear ts_matrix
                    catch
                        
                        timestamps = load(fullfile(fbasename,params.regions{nr},params.eventNames_all{ne},sprintf('Day%s',num2str(days(nf))),flist),variableNames{contains(variableNames,'_ts')});
                        thisVar = variableNames{contains(variableNames,'_ts')|contains(variableNames,'timestamps')};
                        temp.TS = {eval(sprintf('timestamps.%s',thisVar))};
                        clear timestamps
                    end

                    behavEvents = cat(1,behavEvents,temp);
                catch
                    keyboard
                end
            end

        end
            save(fullfile(saveLoc,params.regions{nr},IDs{na},'behavEvents.mat'),'behavEvents','params');

        
    end
end



