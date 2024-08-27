function plotMixedEffects_activeAvoid(temporalKernels,meanFluor,sessIDs,timeBack,timeForward,eventNames,params)



%% Plot each event separately
for nr = 1:numel(params.regions)
    figure('Position',[440 278 1204 420]);

    for ne = 1:numel(eventNames)
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        alpha = linspace(.1,1,numel(sessIDs));
        thisKernel = eval(sprintf('temporalKernels.%s.%s;',params.regions{nr},eventNames{ne}));
        x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
        ax = subplot(2,numel(eventNames),ne); hold on
        plot([0 0], [-.5 1],'--k');
        p  = plot(x,thisKernel,'Color',cmap,'LineWidth',1.5);
        arrayfun(@(x,y) set(p(y),'Color',cat(2,cmap,x)), alpha,sessIDs)
        ylabel('Kernel')
        xlabel('Time from event')
        title(eventNames{ne})
        axis tight

        subplot(2,numel(eventNames),ne+numel(eventNames)); hold on
        for ns = sessIDs
            mu = eval(sprintf('mean(meanFluor.%s.%s(:,:,ns),1,''omitnan'');',params.regions{nr},eventNames{ne}));
            err = eval(sprintf('nansem(meanFluor.%s.%s(:,:,ns),1);',params.regions{nr},eventNames{ne}));
            plot([0 0], [-1 3],'--k');
            plot(x,mu,'Color',cat(2,cmap,alpha(ns)),'LineWidth',1.5);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end
        ylabel('\DeltaF/F');
        xlabel('Time from event')
        axis tight
    end
    
    lh = smallLegend(p,string(sessIDs),gca,'eastoutside');
    lh.Box = 'off';
    lh.Title.String = 'Session';
end


%% Plot all events on the same axis 
for nr = 1:numel(params.regions)
    figure('Position',[440 278 1204 420]);

    % Plot median cross times, cue and shock  
    subplot(2,1,1); hold on    
    title(params.regions{nr})
    plot([0 0], [-1 1],'--k');
    p(3) =  plot([5 5], [-1 1],'Color',params.escapeC);
    plot([params.escapeMed params.escapeMed], [-1 1],'--','Color',[params.escapeC .5]);
    plot([params.avoidMed params.avoidMed], [-1 1],'--','Color',[params.avoidC .5]);
    ylabel('Kernel value');
    xlabel('Time from cue on')

    subplot(2,1,2); hold on
    plot([0 0], [-1 2],'--k');
    plot([5 5], [-1 2],'Color',params.escapeC);
    plot([params.escapeMed params.escapeMed], [-1 2],'--','Color',[params.escapeC .5]);
    plot([params.avoidMed params.avoidMed], [-1 2],'--','Color',[params.avoidC .5]);
    ylabel('\DeltaF/F');
    xlabel('Time from cue on')
   
    for ne = 1:numel(eventNames)

        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        
        % Align x-axis
        thisKernel = eval(sprintf('temporalKernels.%s.%s;',params.regions{nr},eventNames{ne}));
        if contains(eventNames{ne},'Cue')
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
        elseif contains(eventNames{ne}, 'Shock')
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + 5;
        elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Escape')
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + params.escapeMed;
        elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Avoid')
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + params.avoidMed;
        end
        
        subplot(2,1,1);
        p = plot(x,thisKernel,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
        arrayfun(@(x,y) set(p(y),'Color',cat(2,cmap,x)), alpha,sessIDs)
        axis tight

        subplot(2,1,2)
        for ns = sessIDs
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            plot(x,mu,'Color',cat(2,cmap,alpha(ns)),'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end
        axis tight

    end
    
    lh=smallLegend(p,eventNames);
    lh.Title.String = 'Session';
end


%% Plot all events on the same axis
for nr = 1:numel(params.regions)
    figure('Position',[440 41.6667 1204 1.3193e+03]);

    for ns = sessIDs
        % Plot median cross times, cue and shock
        subplot(numel(sessIDs),1,ns); hold on
        if ns == 1
            title(sprintf('%s: session 1',params.regions{nr}))
        else
            title(sprintf('session %s',num2str(ns)))
        end
        plot([0 0], [-1 1],'--k');
        p(3) =  plot([5 5], [-1 1],'Color',params.escapeC);
        plot([params.escapeMed params.escapeMed], [-1 1],'--','Color',[params.escapeC .5]);
        plot([params.avoidMed params.avoidMed], [-1 1],'--','Color',[params.avoidC .5]);
        ylabel('Kernel value');
        if ns == numel(sessIDs)
            xlabel('Time from cue on')
        end

        for ne = 1:numel(eventNames)

            % Plotting parameters based on which event
            if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
                cmap = params.escapeC;
            else
                cmap = params.avoidC;
            end
            if contains(eventNames{ne},'Cross')
                lstyle = '--';
            else
                lstyle = '-';
            end

            % Align x-axis
            thisKernel = eval(sprintf('temporalKernels.%s.%s(ns,:);',params.regions{nr},eventNames{ne}));
            if contains(eventNames{ne},'Cue')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            elseif contains(eventNames{ne}, 'Shock')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + 5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Escape')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + params.escapeMed;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Avoid')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + params.avoidMed;
            end

            p = plot(x,thisKernel,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
            axis tight
        end
    end
end

%% Plot triggered by cue/shock and or cross 
for nr = 1:numel(params.regions)
    figure('Position',[440 278 1204 420]);
    % Plot median cross times, cue and shock
    subplot(2,2,1); hold on
    plot([0 0], [-1 1],'--k');
    l(1)=plot([5 5], [-1 1],'Color',params.escapeC);
    l(2)=plot([params.escapeMed params.escapeMed], [-1 1],'--','Color',[params.escapeC .5]);
    l(3)=plot([params.avoidMed params.avoidMed], [-1 1],'--','Color',[params.avoidC .5]);
    legend(l,{'Shock';'Median escape';'Median avoid'},'AutoUpdate',false)
    
    % Plot cue times
    idx = find(contains(eventNames,'Cue'))';
    for ne = idx
        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        thisKernel = eval(sprintf('temporalKernels.%s.%s;',params.regions{nr},eventNames{ne}));
        x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
        p = plot(x,thisKernel,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
        arrayfun(@(x,y) set(p(y),'Color',cat(2,cmap,x)), alpha,sessIDs)
    end
    ylabel('Kernel value');
    xlabel('Time from cue on')  

    % Plot shock 
    idx = find(contains(eventNames,'Shock'))';
    for ne = idx
        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        thisKernel = eval(sprintf('temporalKernels.%s.%s;',params.regions{nr},eventNames{ne}));
        x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2))+5;
        p = plot(x,thisKernel,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
        arrayfun(@(x,y) set(p(y),'Color',cat(2,cmap,x)), alpha,sessIDs)
    end

    % Plot cross events
    subplot(2,2,2); hold on
    plot([0 0],[-1 1],'--k')
    idx = find(contains(eventNames,'Cross'))';
    for ne = idx
        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        thisKernel = eval(sprintf('temporalKernels.%s.%s;',params.regions{nr},eventNames{ne}));
        x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
        p = plot(x,thisKernel,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
        arrayfun(@(x,y) set(p(y),'Color',cat(2,cmap,x)), alpha,sessIDs)
    end

    ylabel('Kernel value');
    xlabel('Time from cross')


    %%% Plot delta F/F
    subplot(2,2,3); hold on
    plot([0 0], [-1 2],'--k');
    l(1)=plot([5 5], [-1 2],'Color',params.escapeC);
    l(2)=plot([params.escapeMed params.escapeMed], [-1 2],'--','Color',[params.escapeC .5]);
    l(3)=plot([params.avoidMed params.avoidMed], [-1 2],'--','Color',[params.avoidC .5]);

    % Plot cue events
    idx = find(contains(eventNames,'Cue'))';
    for ne = idx
        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        for ns = sessIDs
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p(ne)  = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end
    end
    % Plot shock event
    idx = find(contains(eventNames,'Shock'))';
    for ne = idx
        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        for ns = sessIDs
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2))+5;
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p(3)  = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5);
            patch([x flip(x)],[mu-err flip(mu+err)],p(3).Color,'FaceAlpha',.1,'EdgeColor','none')
        end
    end
    ylabel('\DeltaF/F');
    xlabel('Time from cue on')

    subplot(2,2,4); hold on
    plot([0 0], [-1 2],'--k');
    % Plot cross events
    idx = find(contains(eventNames,'Cross'))';
    for ne = idx
        % Plotting parameters based on which event
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        if contains(eventNames{ne},'Cross')
            lstyle = '--';
        else
            lstyle = '-';
        end
        for ns = sessIDs

            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p(ne)  = plot(x,mu,'--','Color',[cmap alpha(ns)],'LineWidth',1.5);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end
    end
    ylabel('\DeltaF/F');
    xlabel('Time from cross')


end
