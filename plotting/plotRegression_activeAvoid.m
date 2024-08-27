function plotRegression_activeAvoid(temporalKernels,meanFluor,sessIDs,timeBack,timeForward,eventNames,params,saveLoc,fname)



%% Plot each event separately
for nr = 1:numel(params.regions)
    f=figure('Position',[440 278 1204 420]);
    maxY = [];
    minY = [];
    maxY2 = [];
    minY2 = [];
    for ne = 1:numel(eventNames)
        if contains(eventNames{ne},'Escape') || contains(eventNames{ne},'Shock')
            cmap = params.escapeC;
        else
            cmap = params.avoidC;
        end
        alpha = linspace(.1,.95,max(sessIDs));


        subplot(2,numel(eventNames),ne); hold on
     
        for ns = sessIDs
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)); % xaxis for plotting
            mu  = mean(thisKernel,1,'omitnan'); % average kernel value for each session
            err = nansem(thisKernel,1); % sem kernel value for each session
            p(ns)  = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            maxY = cat(1,maxY, max(mu+err));
            minY = cat(1,minY,min(mu-err));
        end
        ylabel('Kernel')
        xlabel('Time from event')
        title(eventNames{ne})
        axis tight

        subplot(2,numel(eventNames),ne+numel(eventNames)); hold on
        maxY2(ne) = max(max(mean(eval(sprintf('meanFluor.%s.%s(:,:,:);',params.regions{nr},eventNames{ne})))));
        minY2(ne) = min(min(mean(eval(sprintf('meanFluor.%s.%s(:,:,:);',params.regions{nr},eventNames{ne})))));
        for ns = sessIDs
            mu = eval(sprintf('mean(meanFluor.%s.%s(:,:,ns),1,''omitnan'');',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(mu,2)); % xaxis for plotting
            err = eval(sprintf('nansem(meanFluor.%s.%s(:,:,ns),1);',params.regions{nr},eventNames{ne}));
            plot(x,mu,'Color',cat(2,cmap,alpha(ns)),'LineWidth',1.5);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            maxY2 = cat(1,maxY2, max(mu+err));
            minY2 = cat(1,minY2,min(mu-err));
        end
        ylabel('\DeltaF/F');
        xlabel('Time from event')
        axis tight
    end

    plotMin(nr) = min(minY);
    plotMax(nr) = max(maxY);
    plotMin2(nr) = min(minY2);
    plotMax2(nr) = max(maxY2);

    for ne = 1:numel(eventNames)
        subplot(2,numel(eventNames),ne); hold on
        plot([0 0], [min(minY) max(maxY)],'--k');
        subplot(2,numel(eventNames),ne+numel(eventNames)); hold on
        plot([0 0], [min(minY2) max(maxY2)], '--k');
    end
    p = p(sessIDs);
    lh = smallLegend(p,string(sessIDs),gca,'eastoutside');
    lh.Box = 'off';
    lh.Title.String = 'Day';
    exportgraphics(f,fullfile(saveLoc,sprintf('%s_%s.pdf',params.regions{nr},fname)))
end


%% Plot all events on the same axis 
for nr = 1:numel(params.regions)
    escapeMed = eval(sprintf('params.escapeMed.%s',params.regions{nr}))+5;
    avoidMed = eval(sprintf('params.avoidMed.%s',params.regions{nr}));

    f=figure('Position',[440 278 1204 420]);
    clear p 
    % Plot median cross times, cue and shock  
    subplot(2,1,1); hold on    
    title(params.regions{nr})
    plot([0 0], [plotMin(nr) plotMax(nr)],'--k');
    p(3) =  plot([5 5],[plotMin(nr) plotMax(nr)],'Color',params.escapeC);
    plot([escapeMed escapeMed], [plotMin(nr) plotMax(nr)],'--','Color',[params.escapeC .5]);
    plot([avoidMed avoidMed], [plotMin(nr) plotMax(nr)],'--','Color',[params.avoidC .5]);
    ylabel('Kernel value');
    xlabel('Time from cue on')

    subplot(2,1,2); hold on
    plot([0 0], [plotMin2(nr) plotMax2(nr)],'--k');
    plot([5 5], [plotMin2(nr) plotMax2(nr)],'Color',params.escapeC);
    plot([escapeMed escapeMed], [plotMin2(nr) plotMax2(nr)],'--','Color',[params.escapeC .5]);
    plot([avoidMed avoidMed], [plotMin2(nr) plotMax2(nr)],'--','Color',[params.avoidC .5]);
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
        for ns = sessIDs
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));

            if contains(eventNames{ne},'Cue')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            elseif contains(eventNames{ne}, 'Shock')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + 5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Escape')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + escapeMed;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Avoid')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + avoidMed;
            end

            subplot(2,1,1);
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p(ne) = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            axis tight

            subplot(2,1,2)
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            plot(x,mu,'Color',cat(2,cmap,alpha(ns)),'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end
        axis tight

    end
    lh=legend(p,eventNames,'Location','northwest','Box','off');
    exportgraphics(f,fullfile(saveLoc,sprintf('%s_%s.pdf',params.regions{nr},fname)),'Append',true)

end


%% Plot all events on the same axis by session
for nr = 1:numel(params.regions)
    f=figure('Position',[433 41.6667 1204 1.3193e+03]);


    for ns = sessIDs
        eval(sprintf('escapeMed = median(params.escapeLatency.%s{ns});',params.regions{nr}))
        eval(sprintf('avoidMed = median(params.avoidLatency.%s{ns});',params.regions{nr}))

        % Plot median cross times, cue and shock
        subplot(max(params.sessIDs),2,ns+(ns-1)); hold on
        title(sprintf('%s : Day %s',params.regions{nr},num2str(ns)))
        plot([0 0], [plotMin(nr) plotMax(nr)],'--k');
        p(3) =  plot([5 5],[plotMin(nr) plotMax(nr)],'Color',params.escapeC);
        plot([escapeMed+5 escapeMed+5], [plotMin(nr) plotMax(nr)],'--','Color',[params.escapeC .5]);
        plot([avoidMed avoidMed], [plotMin(nr) plotMax(nr)],'--','Color',[params.avoidC .5]);
        ylabel('Kernel value');
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
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu  = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            if contains(eventNames{ne},'Cue')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            elseif contains(eventNames{ne}, 'Shock')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + 5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Escape')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + escapeMed +5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Avoid')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + avoidMed;
            end

            p = plot(x,mu,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            axis tight
        end

        % Plot fluorescence
        % Plot median cross times, cue and shock
        subplot(max(params.sessIDs),2,ns*2); hold on
        title(sprintf('%s : Day %s',params.regions{nr},num2str(ns)))
        plot([0 0], [plotMin2(nr) plotMax2(nr)],'--k');
        p(3) =  plot([5 5],[plotMin2(nr) plotMax2(nr)],'Color',params.escapeC);
        plot([escapeMed+5 escapeMed+5], [plotMin2(nr) plotMax2(nr)],'--','Color',[params.escapeC .5]);
        plot([avoidMed avoidMed], [plotMin2(nr) plotMax2(nr)],'--','Color',[params.avoidC .5]);
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
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu  = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            if contains(eventNames{ne},'Cue')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            elseif contains(eventNames{ne}, 'Shock')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + 5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Escape')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + escapeMed + 5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Avoid')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + avoidMed;
            end

            p = plot(x,mu,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            axis tight
        end

    end
    exportgraphics(f,fullfile(saveLoc,sprintf('%s_%s.pdf',params.regions{nr},fname)),'Append',true)
end

%% Plot triggered by cue/shock and or cross 
for nr = 1:numel(params.regions)
    escapeMed = eval(sprintf('params.escapeMed.%s',params.regions{nr}));
    avoidMed = eval(sprintf('params.avoidMed.%s',params.regions{nr}));

 
    f=figure('Position',[440 278 1204 420]);
    % Plot median cross times, cue and shock
    subplot(2,2,1); hold on
    plot([0 0], [plotMin(nr) plotMax(nr)],'--k');
    l(1)=plot([5 5],[plotMin(nr) plotMax(nr)],'Color',params.escapeC);
    l(2)=plot([escapeMed+5 escapeMed+5], [plotMin(nr) plotMax(nr)],'--','Color',[params.escapeC .5]);
    l(3)=plot([avoidMed avoidMed], [plotMin(nr) plotMax(nr)],'--','Color',[params.avoidC .5]);
    legend(l,{'Shock';'Median escape';'Median avoid'},'AutoUpdate',false,'Box','off')
    subplot(2,2,2); hold on
    plot([0 0], [plotMin(nr) plotMax(nr)],'--k')
    % Plot cue times
    for ns = sessIDs
        subplot(2,2,1); hold on
        idx = contains(eventNames,'Cue');
        for ne = find(idx)'
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
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            p = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end

        ylabel('Kernel value');
        xlabel('Time from cue on')
        axis tight

        % Plot shock
        idx = find(contains(eventNames,'Shock'));
        for ne = idx'
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
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2))+5;
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end

        % Plot cross events
        subplot(2,2,2); hold on
        idx = find(contains(eventNames,'Cross'));
        for ne = idx'
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
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end

        ylabel('Kernel value');
        xlabel('Time from cross')
        axis tight
    end

    %%% Plot delta F/F
    subplot(2,2,3); hold on
    plot([0 0], [plotMin2(nr) plotMax2(nr)],'--k');
    l(1)=plot([5 5],  [plotMin2(nr) plotMax2(nr)],'Color',params.escapeC);
    l(2)=plot([escapeMed+5 escapeMed+5],  [plotMin2(nr) plotMax2(nr)],'--','Color',[params.escapeC .5]);
    l(3)=plot([avoidMed avoidMed], [plotMin2(nr) plotMax2(nr)],'--','Color',[params.avoidC .5]);
    subplot(2,2,4); hold on
    plot([0 0], [plotMin2(nr) plotMax2(nr)],'--k')
    % Plot cue times
    for ns = sessIDs
        subplot(2,2,3); hold on
        idx = contains(eventNames,'Cue');
        for ne = find(idx)'
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
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            p = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end

        ylabel('\DeltaF/F');
        xlabel('Time from cross')
        axis tight

        % Plot shock
        idx = find(contains(eventNames,'Shock'));
        for ne = idx'
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
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2))+5;
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end

        % Plot cross events
        subplot(2,2,4); hold on
        idx = find(contains(eventNames,'Cross'));
        for ne = idx'
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
            thisKernel = eval(sprintf('meanFluor.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            mu = mean(thisKernel,1,'omitnan');
            err = nansem(thisKernel,1);
            p = plot(x,mu,'Color',[cmap alpha(ns)],'LineWidth',1.5,'LineStyle',lstyle);
            patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
        end

        ylabel('\DeltaF/F');
        xlabel('Time from cross')
        axis tight
    end
    exportgraphics(f,fullfile(saveLoc,sprintf('%s_%s.pdf',params.regions{nr},fname)),'Append',true)


end
