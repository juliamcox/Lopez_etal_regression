function  plotGoodnessOfFit_activeAvoid(eventNames,sessIDs,modelStats,saveLoc,fname,params,modelStats_shuff)

al = .05/numel(eventNames); 


if isfield(modelStats{1,1},'Fstat') && params.MEFlag == 0

    %% Plot proportion of mice with significant kernels (if fit by animal)
    if size(modelStats,1)>1
        if params.numShuff > 0
            f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
            t=tiledlayout(numel(eventNames),1);
            for ne = 1:numel(eventNames)
                clear thisP
                for na = 1:size(modelStats,1)
                    for ns = 1:numel(sessIDs)
                        if ~isnan(modelStats{na,ns}.corr.full)
                            eval(sprintf('thisP(na,ns) = modelStats{na,ns}.pvals.Fstat.%s;',eventNames{ne}));
                        else
                            thisP(na,ns) = NaN;
                        end
                    end
                end
                nexttile
                hold on
                mu  = mean(thisP<al,1,'omitnan');
                bar(mu,'FaceColor','k','EdgeColor','none','FaceAlpha',.5)
                title(eventNames{ne})
                xlabel(t,'Day')
                ylabel(t,'P(significant)')
                set(gca,'XTick',sessIDs,'YLim',[0 1])
            end
            exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

            % no refit
            f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
            t=tiledlayout(numel(eventNames),1);
            for ne = 1:numel(eventNames)
                clear thisP
                for na = 1:size(modelStats,1)
                    for ns = 1:numel(sessIDs)
                        if ~isnan(modelStats{na,ns}.corr.full)
                            eval(sprintf('thisP(na,ns) = modelStats{na,ns}.pvals.Fstat_noRefit.%s;',eventNames{ne}));
                        else
                            thisP(na,ns) = NaN;
                        end
                    end
                end
                nexttile
                hold on
                mu  = mean(thisP<al,1,'omitnan');
                bar(mu,'FaceColor','k','EdgeColor','none','FaceAlpha',.5)
                title(eventNames{ne})
                xlabel(t,'Day')
                ylabel(t,'P(significant w/o refitting)')
                set(gca,'XTick',sessIDs,'YLim',[0 1])
            end
            exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)
        end
    end

    %% Plot F-statistic for full vs. reduced models

    f=figure('Position',[653 453.6667 1.0913e+03 563.3333]);
    t=tiledlayout(1,numel(params.sessIDs));
    plotmin = 0;
    for ns = sessIDs
        p(ns)=nexttile;
        hold on
        clear thisF thisP
        thisShuff = cell(size(eventNames));
        for ne = 1:numel(eventNames)
            for na = 1:size(modelStats,1)
                if isnan(modelStats{na,ns}.corr.full)
                    thisF(na,ne) = NaN;
                else
                    try
                        eval(sprintf('thisF(na,ne) = modelStats{na,ns}.Fstat.%s;',eventNames{ne}));
                    catch
                        thisF(na,ne) = NaN;
                    end
                end
            end

            % If fit by animal, perform a 1 sample t-test for each event versus 0 or chance (as determine w/ shuffle)

            if params.numShuff > 0 % if ran shuffles
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    if ~isnan(modelStats{na,ns}.corr.full)
                        for nss = 1:params.numShuff
                            tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.Fstat.%s;',eventNames{ne}));
                        end
                        thisShuff{ne} = cat(2,thisShuff{ne},tempShuff);
                    end
                end
            end
            if size(modelStats,1)>1
                if params.numShuff >0
                    [~,thisP(ne)] = ttest(thisF(:,ne),mean(thisShuff{ne},2,'omitnan'));
                else
                    [~,thisP(ne)] = ttest(thisF(:,ne));
                end
            else
                thisP(ne) = false(1,1);
            end
        end
        mu = mean(thisF,1,'omitnan');
        err= nansem(thisF,1);
        x = 1:numel(mu);
        scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k');
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color','k','LineWidth',1)
        ymax(ns) = max(mu+err);
        
        scatter(x(thisP<al),mu(thisP<al),20,'MarkerFaceColor','k','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
        box off
        title(sprintf('Day %s',num2str(ns)))
        set(gca,'XTick',x)
        set(gca, 'XTickLabel',eventNames,'XTickLabelRotation',50)
        % If shuffles, plot mean + standard deviation of the shuffle
        
        if params.numShuff > 0
            for ne = 1:numel(eventNames)
                yup(ne)  = mean(thisShuff{ne},2,'omitnan') + std(thisShuff{ne},[],'omitnan');
                ylow(ne) = mean(thisShuff{ne},2,'omitnan') - std(thisShuff{ne},[],'omitnan');
                plot([x(ne)-.25 x(ne)+.25],[yup(ne) yup(ne)],'k')
                plot([x(ne)-.25 x(ne)+.25],[ylow(ne) ylow(ne)],'k')
            end
            plotmin = cat(2,plotmin,ylow);
        else 
            plotmin = 0; 
        end
       
        
    end
    xlabel(t,'Removed event')
    ylabel(t,'F-statistic (Full vs. reduced)')
    for ns = 1:numel(sessIDs)
        set(p(ns),'YLim',[min(plotmin)-.5 max(ymax)+.5],'XLim',[0 numel(eventNames)+1])
    end
    exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

    % no refit
    f=figure('Position',[653 453.6667 1.0913e+03 563.3333]);
    t=tiledlayout(1,numel(params.sessIDs));
    plotmin = 0;
    clear ymax
    for ns = sessIDs
        p(ns)=nexttile;
        hold on
        clear thisF thisP
        thisShuff = cell(size(eventNames));
        for ne = 1:numel(eventNames)
            for na = 1:size(modelStats,1)
                if isnan(modelStats{na,ns}.corr.full)
                    thisF(na,ne) = NaN;
                else
                    try
                        eval(sprintf('thisF(na,ne) = modelStats{na,ns}.Fstat_noRefit.%s;',eventNames{ne}));
                    catch
                        thisF(na,ne) = NaN;
                    end
                end
            end

            % If fit by animal, perform a 1 sample t-test for each event versus 0 or chance (as determine w/ shuffle)
              if params.numShuff > 0 % if ran shuffles
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    if ~isnan(modelStats{na,ns}.corr.full)
                        for nss = 1:params.numShuff
                            tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.Fstat_noRefit.%s;',eventNames{ne}));
                        end
                        thisShuff{ne} = cat(2,thisShuff{ne},tempShuff);
                    end
                end
            end
            if size(modelStats,1)>1
                if params.numShuff >0
                    [~,thisP(ne)] = ttest(thisF(:,ne),mean(thisShuff{ne},2,'omitnan'));
                else
                    [~,thisP(ne)] = ttest(thisF(:,ne));
                end
            else
                thisP(ne) = false(1,1);
            end
        end
        mu = mean(thisF,1,'omitnan');
        err= nansem(thisF,1);
        x = 1:numel(mu);
        scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k');
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color','k','LineWidth',1)
        ymax(ns) = max(mu+err);
        scatter(x(thisP<al),mu(thisP<al),20,'MarkerFaceColor','k','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
        box off
        title(sprintf('Day %s',num2str(ns)))
        set(gca,'XTick',x)
        set(gca, 'XTickLabel',eventNames,'XTickLabelRotation',50)

        % If shuffles, plot mean + standard deviation of the shuffle
        if params.numShuff > 0
            for ne = 1:numel(eventNames)
                yup(ne)  = mean(thisShuff{ne},2,'omitnan') + std(thisShuff{ne},[],'omitnan');
                ylow(ne) = mean(thisShuff{ne},2,'omitnan') - std(thisShuff{ne},[],'omitnan');
                plot([x(ne)-.25 x(ne)+.25],[yup(ne) yup(ne)],'k')
                plot([x(ne)-.25 x(ne)+.25],[ylow(ne) ylow(ne)],'k')
            end
             plotmin = cat(2,plotmin,ylow);
        else
            plotmin = 0;
        end
       
    end
    xlabel(t,'Removed event (without refitting)')
    ylabel(t,'F-statistic (Full vs. reduced)')
    for ns = 1:numel(sessIDs)
        set(p(ns),'YLim',[min(plotmin)-.5 max(ymax)+.5],'XLim',[0 numel(eventNames)+1])
    end
    exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

    % Plot by event 
    
    f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
    t=tiledlayout(numel(eventNames),1);
    plotmin = 0;
    clear ymax
    for ne = 1:numel(eventNames)
        clear thisF thisP
        thisShuff = cell(size(sessIDs));
        p(ne)=nexttile;
        hold on
        for ns=sessIDs
            for na = 1:size(modelStats,1)
                if isnan(modelStats{na,ns}.corr.full)
                    thisF(na,ns) = NaN;
                else
                    eval(sprintf('thisF(na,ns) = modelStats{na,ns}.Fstat.%s;',eventNames{ne}));
                end
            end
            % If fit by animal, perform a 1 sample t-test for each event versus 0 or chance (as determine w/ shuffle)
            if params.numShuff > 0 % if ran shuffles
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    if ~isnan(modelStats{na,ns}.corr.full)
                        for nss = 1:params.numShuff
                            tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.Fstat.%s;',eventNames{ne}));
                        end
                        thisShuff{ns} = cat(2,thisShuff{ns},tempShuff);
                    end
                end
            end
            if size(modelStats,1)>1
                if params.numShuff >0
                    [~,thisP(ns)] = ttest(thisF(:,ns),mean(thisShuff{ns},2,'omitnan'));
                else
                    [~,thisP(ns)] = ttest(thisF(:,ns));
                end
            else
                thisP(ns) = false(1,1);
            end
        end

        mu = mean(thisF,1,'omitnan');
        err= nansem(thisF,1);
        x = 1:numel(mu);
        scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k');
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color','k','LineWidth',1)
        ymax(ne) = max(mu+err);
        scatter(x(thisP<al),mu(thisP<al),20,'MarkerFaceColor','k','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
        box off


        % If fit by animal, 1 way anova by day
        if size(modelStats,1)>1
            aov = anova(thisF);
            aov = stats(aov);
            title(sprintf('%s: 1-way ANOVA p = %.4f',eventNames{ne},(aov.pValue)))
        else
            title(eventNames{ne})
        end

        set(gca,'XTick',x)

        % If shuffles, plot mean + standard deviation of the shuffle
        if params.numShuff > 0
            for ns = 1:numel(sessIDs)
                yup(ns)  = mean(thisShuff{ns},2,'omitnan') + std(thisShuff{ns},[],'omitnan');
                ylow(ns) = mean(thisShuff{ns},2,'omitnan') - std(thisShuff{ns},[],'omitnan');
                plot([x(ns)-.25 x(ns)+.25],[yup(ns) yup(ns)],'k')
                plot([x(ns)-.25 x(ns)+.25],[ylow(ns) ylow(ns)],'k')
                plotmin = cat(2,plotmin,ylow);
            end
        else
            plotmin = 0;
        end
    end
    xlabel(t,'Day')
    ylabel(t,'F-statistic')
    for ne = 1:numel(eventNames)
        set(p(ne),'YLim',[min(plotmin)-.5 max(ymax)+.5],'XLim',[0 numel(sessIDs)+1])
    end
    exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

    % Plot by event no refit
    f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
    t=tiledlayout(numel(eventNames),1);
    plotmin = 0;
    clear ymax
    for ne = 1:numel(eventNames)
        clear thisF thisP
        thisShuff = cell(size(sessIDs));
        p(ne)=nexttile;
        hold on
        for ns=sessIDs
            for na = 1:size(modelStats,1)
                if isnan(modelStats{na,ns}.corr.full)
                    thisF(na,ns) = NaN;
                else
                    eval(sprintf('thisF(na,ns) = modelStats{na,ns}.Fstat_noRefit.%s;',eventNames{ne}));
                end
            end
            % If fit by animal, perform a 1 sample t-test for each event versus 0 or chance (as determined w/ shuffle)
            if params.numShuff > 0 % if ran shuffles
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    if ~isnan(modelStats{na,ns}.corr.full)
                        for nss = 1:params.numShuff
                            tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.Fstat_noRefit.%s;',eventNames{ne}));
                        end
                        thisShuff{ns} = cat(2,thisShuff{ns},tempShuff);
                    end
                end
            end
            if size(modelStats,1)>1
                if params.numShuff >0
                    [~,thisP(ns)] = ttest(thisF(:,ns),mean(thisShuff{ns},2,'omitnan'));
                else
                    [~,thisP(ns)] = ttest(thisF(:,ns));
                end
            else
                thisP(ns) = false(1,1);
            end
        end

        mu = mean(thisF,1,'omitnan');
        err= nansem(thisF,1);
        x = 1:numel(mu);
        scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k');
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color','k','LineWidth',1)
        ymax(ne) = max(mu+err);
        scatter(x(thisP<al),mu(thisP<al),20,'MarkerFaceColor','k','MarkerFaceAlpha',.5,'MarkerEdgeColor','none')
        box off
        if size(modelStats,1)>1
            aov = anova(thisF);
            aov = stats(aov);
            title(sprintf('%s: 1-way ANOVA p = %.4f',eventNames{ne},(aov.pValue)))
        else
            title(eventNames{ne})
        end       
        set(gca,'XTick',x)

        % If shuffles, plot mean + standard deviation of the shuffle
        if params.numShuff > 0
            for ns = 1:numel(sessIDs)
                yup(ns)  = mean(thisShuff{ns},2,'omitnan') + std(thisShuff{ns},[],'omitnan');
                ylow(ns) = mean(thisShuff{ns},2,'omitnan') - std(thisShuff{ns},[],'omitnan');
                plot([x(ns)-.25 x(ns)+.25],[yup(ns) yup(ns)],'k')
                plot([x(ns)-.25 x(ns)+.25],[ylow(ns) ylow(ns)],'k')
            plotmin = cat(2,plotmin,ylow);
            end
        else
            plotmin = 0;
        end
    end
    xlabel(t,'Day')
    ylabel(t,'F-statistic (no refitting)')
    for ne = 1:numel(eventNames)
        set(p(ne),'YLim',[min(plotmin)-.5 max(ymax)+.5],'XLim',[0 numel(sessIDs)+1])
    end
    exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

end

%% Plot correlation coefficient for full and reduced models


f=figure('Position',[653 453.6667 1.0913e+03 563.3333]);
t=tiledlayout(1,numel(params.sessIDs));
plotmin = 0;
clear ymax
for ns = sessIDs
    thisShuff =cell(numel(eventNames)+1);
    nexttile
    hold on
    clear thisF thisP
    for na = 1:size(modelStats,1)
        if isnan(modelStats{na,ns}.corr.full)
            thisF(na,1) = NaN;
        else
            thisF(na,1) = modelStats{na,ns}.corr.full;
        end
    end
    for ne = 1:numel(eventNames)
        for na = 1:size(modelStats,1)
            if isnan(modelStats{na,ns}.corr.full)
                thisF(na,ne+1) = NaN;
            else
                eval(sprintf('thisF(na,ne+1) = modelStats{na,ns}.corr.%s;',eventNames{ne}));
            end
        end
    end


    % If fit by animal, perform a 1 sample t-test for each event versus 0 or chance (as determine w/ shuffle)
    if size(modelStats,1)>1
        if params.numShuff > 0 % if ran shuffles
            for na = 1:size(modelStats,1)
                clear tempShuff
                if ~isnan(modelStats{na,ns}.corr.full)
                    for nss = 1:params.numShuff
                        tempShuff(nss) = modelStats_shuff{na,ns}{nss}.corr.full;
                    end
                    thisShuff{1} = cat(2,thisShuff{1},tempShuff);
                end
            end
            [~,thisP(1)] = ttest(thisF(:,1),mean(thisShuff{1},2,'omitnan'));
            for ne = 1:numel(eventNames)
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    if ~isnan(modelStats{na,ns}.corr.full)
                        for nss = 1:params.numShuff
                            tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.corr.%s;',eventNames{ne}));
                        end
                        thisShuff{ne+1} = cat(2,thisShuff{ne+1},tempShuff);
                    end
                end
                [~,thisP(ne+1)] = ttest(thisF(:,ne),mean(thisShuff{ne},2,'omitnan'));
            end
        else
            [~,thisP] = ttest(thisF);
        end
    else
        if params.numShuff>0
        for na = 1:size(modelStats,1)
            clear tempShuff
            for nss = 1:params.numShuff
                tempShuff(nss) = modelStats_shuff{na,ns}{nss}.corr.full;
            end
            thisShuff{1} = cat(2,thisShuff{1},tempShuff);
        end
        for ne = 1:numel(eventNames)
            for na = 1:size(modelStats,1)
                clear tempShuff
                for nss = 1:params.numShuff
                    tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.corr.%s;',eventNames{ne}));
                end
                thisShuff{ne+1} = cat(2,thisShuff{ne+1},tempShuff);
            end
        end
        end
        thisP  = nan(numel(eventNames)+1);
        
    end

    mu = mean(thisF,1,'omitnan');
    err= nansem(thisF,1);
    x = 1:numel(mu);
    scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
    errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
    ymax(ns) = max(mu+err);
    sigFlag = thisP<al;
    scatter(x(sigFlag),mu(sigFlag),20,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.5)

    if params.numShuff > 0
        % If shuffles, plot mean + standard deviation of the shuffle
        for ne = 1:numel(eventNames)+1
            yup(ne)  = mean(thisShuff{ne},2,'omitnan') + std(thisShuff{ne},[],'omitnan');
            ylow(ne) = mean(thisShuff{ne},2,'omitnan') - std(thisShuff{ne},[],'omitnan');
            plot([x(ne)-.25 x(ne)+.25],[yup(ne) yup(ne)],'k')
            plot([x(ne)-.25 x(ne)+.25],[ylow(ne) ylow(ne)],'k')
        end
        plotmin = cat(2,plotmin,ylow);
    end
    title(sprintf('Day %s',num2str(ns)))
    set(gca,'XTick', x, 'XTickLabel',cat(1,{'Full'},cellfun(@(x) cat(2,'- ', x),eventNames,'UniformOutput',false)),'XTickLabelRotation',90,'XLim',[0 max(x)+1])
end
p = t.Children;
for ns = sessIDs
    set(p(ns),'YLim',[-max(ymax)+.3 max(ymax)+.1])
end

ylabel(t,'Corr(real,estimated)')
xlabel(t,'Model')
exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

% no refit
f=figure('Position',[653 453.6667 1.0913e+03 563.3333]);
t=tiledlayout(1,numel(params.sessIDs));
plotmin = 0;
clear ymax
for ns = sessIDs
    thisShuff =cell(numel(eventNames)+1);
    nexttile
    hold on
    clear thisF thisP
    for na = 1:size(modelStats,1)
        if isnan(modelStats{na,ns}.corr.full)
            thisF(na,1) = NaN;
        else
            thisF(na,1) = modelStats{na,ns}.corr.full;
        end
    end
    for ne = 1:numel(eventNames)
        for na = 1:size(modelStats,1)
            if isnan(modelStats{na,ns}.corr.full)
                thisF(na,ne+1) = NaN;
            else
                eval(sprintf('thisF(na,ne+1) = modelStats{na,ns}.corr_noRefit.%s;',eventNames{ne}));
            end
        end
    end


    % If fit by animal, perform a 1 sample t-test for each event versus 0 or chance (as determine w/ shuffle)
    if size(modelStats,1)>1
        if params.numShuff > 0 % if ran shuffles
            for na = 1:size(modelStats,1)
                clear tempShuff
                if ~isnan(modelStats{na,ns}.corr.full)
                    for nss = 1:params.numShuff
                        tempShuff(nss) = modelStats_shuff{na,ns}{nss}.corr.full;
                    end
                    thisShuff{1} = cat(2,thisShuff{1},tempShuff);
                end
            end
            [~,thisP(1)] = ttest(thisF(:,1),mean(thisShuff{1},2,'omitnan'));
            for ne = 1:numel(eventNames)
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    if ~isnan(modelStats{na,ns}.corr.full)
                        for nss = 1:params.numShuff
                            tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.corr_noRefit.%s;',eventNames{ne}));
                        end
                        thisShuff{ne+1} = cat(2,thisShuff{ne+1},tempShuff);
                    end
                end
                [~,thisP(ne+1)] = ttest(thisF(:,ne),mean(thisShuff{ne+1},2,'omitnan'));
            end
        else
            [~,thisP] = ttest(thisF);
        end
    else
        if params.numShuff>0
            for na = 1:size(modelStats,1)
                clear tempShuff
                for nss = 1:params.numShuff
                    tempShuff(nss) = modelStats_shuff{na,ns}{nss}.corr.full;
                end
                thisShuff{1} = cat(2,thisShuff{1},tempShuff);
            end
            for ne = 1:numel(eventNames)
                for na = 1:size(modelStats,1)
                    clear tempShuff
                    for nss = 1:params.numShuff
                        tempShuff(nss) = eval(sprintf('modelStats_shuff{na,ns}{nss}.corr_noRefit.%s;',eventNames{ne}));
                    end
                    thisShuff{ne+1} = cat(2,thisShuff{ne+1},tempShuff);
                end
            end
        end
        thisP = nan(numel(eventNames)+1);
    end

    mu = mean(thisF,1,'omitnan');
    err= nansem(thisF,1);
    x = 1:numel(mu);
    scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
    errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
    ymax(ns) = max(mu+err);
    sigFlag = thisP<al;
    scatter(x(sigFlag),mu(sigFlag),20,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.5)

    if params.numShuff > 0
        % If shuffles, plot mean + standard deviation of the shuffle
        for ne = 1:numel(eventNames)+1
            yup(ne)  = mean(thisShuff{ne},2,'omitnan') + std(thisShuff{ne},[],'omitnan');
            ylow(ne) = mean(thisShuff{ne},2,'omitnan') - std(thisShuff{ne},[],'omitnan');
            plot([x(ne)-.25 x(ne)+.25],[yup(ne) yup(ne)],'k')
            plot([x(ne)-.25 x(ne)+.25],[ylow(ne) ylow(ne)],'k')
        end
        plotmin = cat(2,plotmin,ylow);
    end
    title(sprintf('Day %s',num2str(ns)))
    set(gca,'XTick', x, 'XTickLabel',cat(1,{'Full'},cellfun(@(x) cat(2,'- ', x),eventNames,'UniformOutput',false)),'XTickLabelRotation',90,'XLim',[0 max(x)+1])
end
p = t.Children;
for ns = sessIDs
    set(p(ns),'YLim',[-max(ymax)+.3 max(ymax)+.1])
end

ylabel(t,'Corr(real,estimated - no refitting)')
xlabel(t,'Model')
exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

%% Plot correlation coefficient for full and reduced models by event

f=figure('Position', [639.6667 843 701.3333 476.0333]);
t=tiledlayout(2,1);
clear thisF 
% Plot correlation coefficient for full model across days
for na = 1:size(modelStats,1)
    for ns = sessIDs
        if isnan(modelStats{na,ns}.corr.full)
            thisF(na,ns) = NaN;
        else
            thisF(na,ns) = modelStats{na,ns}.corr.full;
        end
    end
end
nexttile;
hold on
mu = mean(thisF,1,'omitnan');
err= nansem(thisF,1);
x = 1:numel(mu);
scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
ylabel('Corr(real,estimated)')
xlabel('Day')
box off
ylim = get(gca,'YLim');
set(gca,'XLim', [0 numel(sessIDs)+1],'YLim',[0 ylim(2)],'XTick',sessIDs);
title('Full model')

% 1-way ANOVA by session
aov = anova(thisF,'FactorNames',"Day");
thisP = stats(aov);
%columname = this
pPlot = table2cell(thisP);
pPlot = cat(2,thisP.Properties.RowNames,pPlot);
pPlot = cat(1,cat(2,{' '},thisP.Properties.VariableNames),pPlot);
nexttile
hold on 
g = gca;
g.XAxis.Color = 'none';
g.YAxis.Color = 'none';
g.Color = 'none';
for nr = 1:size(pPlot,1)
    text((ones(size(pPlot,2),1)+1.5.*[1:size(pPlot,2)]'),ones(size(pPlot,2),1)-.25*nr,pPlot(nr,:),'FontSize',10)
end
g.XLim = [3.5 10];




exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)




%% Comparison between correlation coefficient for full and reduced models by session for each event 

f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
t=tiledlayout(numel(eventNames),1);

for ne = 1:numel(eventNames)
    clear thisF ymax thisEvent
    % Plot correlation coefficient for full model across days
    for na = 1:size(modelStats,1)
        for ns = sessIDs
            if isnan(modelStats{na,ns}.corr.full)
                thisF(na,ns) = NaN;
                thisEvent(na,ns) = NaN;
            else
                thisF(na,ns) = modelStats{na,ns}.corr.full;
                thisEvent(na,ns) = eval(sprintf('modelStats{na,ns}.corr.%s',eventNames{ne}));
            end
        end
    end
    [~,thisP] = ttest2(thisF,thisEvent);
    sigFlag = thisP<(.05/numel(sessIDs));
    
    nexttile;
    hold on
    mu_full  = mean(thisF,1,'omitnan');
    err_full = nansem(thisF,1);
    mu      = mean(thisEvent,1,'omitnan');
    err     = nansem(thisEvent,1);
    x       = (1:numel(mu))+.1;
    x_full  = (1:numel(mu))-.1;
    %scatter(x,mu,20,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
    if size(modelStats,1)>1
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
        errorbar(x_full,mu_full,err_full,'LineStyle','none','CapSize',0,'Color',[72/255 209/255 204/255 .4],'LineWidth',1);
    else
        scatter(x,mu,10,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
        scatter(x_full,mu_full,10,'MarkerFaceColor','w','MarkerEdgeColor',[72/255 209/255 204/255],'MarkerFaceAlpha',.5)
    end
    
    
    text(x(sigFlag)-.1,(err_full(sigFlag)+mu_full(sigFlag)+.1),'*')
    box off
    ylim = get(gca,'YLim');
    set(gca,'XLim', [0 numel(sessIDs)+1],'YLim',[0 ylim(2)]);
    title(eventNames{ne})
    ymax(ne) = max(err_full+mu_full)+.15;
    
end

p = t.Children;
for ne = 1:numel(p)
    p(ne).YLim = [-max(ymax)+.5 max(ymax)+.1];
end
ylabel(t,'Corr(real,estimated)')
xlabel(t,'Day')
exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
t=tiledlayout(numel(eventNames),1);

for ne = 1:numel(eventNames)
    clear thisF ymax thisEvent
    % Plot correlation coefficient for full model across days
    for na = 1:size(modelStats,1)
        for ns = sessIDs
            if isnan(modelStats{na,ns}.corr.full)
                thisF(na,ns) = NaN;
                thisEvent(na,ns) = NaN;
            else
                thisF(na,ns) = modelStats{na,ns}.corr.full;
                thisEvent(na,ns) = eval(sprintf('modelStats{na,ns}.corr_noRefit.%s',eventNames{ne}));
            end
        end
    end
    if size(modelStats,1)>1
    [~,thisP] = ttest2(thisF,thisEvent);
    sigFlag = thisP<(.05/numel(sessIDs));
    else
        thisP = false(size(sessIDs));
    end
    nexttile;
    hold on
    mu_full  = mean(thisF,1,'omitnan');
    err_full = nansem(thisF,1);
    mu      = mean(thisEvent,1,'omitnan');
    err     = nansem(thisEvent,1);
    x       = (1:numel(mu))+.1;
    x_full  = (1:numel(mu))-.1;
    if size(modelStats,1)>1
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
        errorbar(x_full,mu_full,err_full,'LineStyle','none','CapSize',0,'Color',[72/255 209/255 204/255 .4],'LineWidth',1);
    else
        scatter(x,mu,10,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
        scatter(x_full,mu_full,10,'MarkerFaceColor','w','MarkerEdgeColor',[72/255 209/255 204/255],'MarkerFaceAlpha',.5)
    end
    text(x(sigFlag)-.1,(err_full(sigFlag)+mu_full(sigFlag)+.1),'*')
    box off
    ylim = get(gca,'YLim');
    set(gca,'XLim', [0 numel(sessIDs)+1],'YLim',[0 ylim(2)]);
    title(eventNames{ne})
    ymax(ne) = max(err_full+mu_full)+.15;
    
end

p = t.Children;
for ne = 1:numel(p)
    p(ne).YLim = [-max(ymax)+.5 max(ymax)+.1];
end
ylabel(t,'Corr(real,estimated) - no refitting')
xlabel(t,'Day')
exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

%% Plot diff between full and reduced model across days 

f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
t=tiledlayout(numel(eventNames),1);

for ne = 1:numel(eventNames)
    clear thisF ymax 
    % Plot correlation coefficient for full model across days
    for na = 1:size(modelStats,1)
        for ns = sessIDs
            if isnan(modelStats{na,ns}.corr.full)
                thisF(na,ns) = NaN;
            else
                thisF(na,ns) = modelStats{na,ns}.corr.full-eval(sprintf('modelStats{na,ns}.corr.%s',eventNames{ne}));
            end
        end
    end
    
    nexttile;
    hold on
   
    mu      = mean(thisF,1,'omitnan');
    err     = nansem(thisF,1);
    x       = (1:numel(mu));
    if size(modelStats,1)>1
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
    else
        scatter(x,mu,10,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
    end
    box off
    ylim = get(gca,'YLim');
    try
        set(gca,'XLim', [0 numel(sessIDs)+1],'YLim',[0 ylim(2)]);
    catch
        set(gca,'XLim', [0 numel(sessIDs)+1],'YLim',[ylim(1) 0]);
    end

    title(eventNames{ne})
    ymax(ne) = max(err+mu);
    plot([0 numel(sessIDs)+1],[0 0], '--','Color',[0 0 0 .4])
    if size(modelStats,1)>1
        aov = anova(thisF);
        aov = stats(aov);
        title(sprintf('%s: 1-way ANOVA p = %.4f',eventNames{ne},(aov.pValue)))
    else
        title(eventNames{ne})
    end
    
end

% p = t.Children;
% for ne = 1:numel(p)
%     p(ne).YLim = [-max(ymax) max(ymax)+.01];
% end
ylabel(t,'Corr coefficient: Full-Reduced')
xlabel(t,'Day')
exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

% no refit 
f=figure('Position', [639.6667 268.3333 467.3333 1.0507e+03]);
t=tiledlayout(numel(eventNames),1);

for ne = 1:numel(eventNames)
    clear thisF ymax 
    % Plot correlation coefficient for full model across days
    for na = 1:size(modelStats,1)
        for ns = sessIDs
            if isnan(modelStats{na,ns}.corr.full)
                thisF(na,ns) = NaN;
            else
                thisF(na,ns) = modelStats{na,ns}.corr.full-eval(sprintf('modelStats{na,ns}.corr_noRefit.%s',eventNames{ne}));
            end
        end
    end
    
    nexttile;
    hold on
   
    mu      = mean(thisF,1,'omitnan');
    err     = nansem(thisF,1);
    x       = (1:numel(mu));
    if size(modelStats,1)>1
        errorbar(x,mu,err,'LineStyle','none','CapSize',0,'Color',[0 0 0 .4],'LineWidth',1);
    else
        scatter(x,mu,10,'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
    end
    box off
    ylim = get(gca,'YLim');
    set(gca,'XLim', [0 numel(sessIDs)+1],'YLim',[0 ylim(2)]);
    title(eventNames{ne})
    ymax(ne) = max(err+mu);
    plot([0 numel(sessIDs)+1],[0 0], '--','Color',[0 0 0 .4])
    if size(modelStats,1)>1
        aov = anova(thisF);
        aov = stats(aov);
        title(sprintf('%s: 1-way ANOVA p = %.4f',eventNames{ne},(aov.pValue)))
    else
        title(eventNames{ne})
    end
    
end

% p = t.Children;
% for ne = 1:numel(p)
%     p(ne).YLim = [-max(ymax) max(ymax)+.01];
% end
ylabel(t,'Corr coefficient: Full-Reduced, no refitting')
xlabel(t,'Day')
exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)
