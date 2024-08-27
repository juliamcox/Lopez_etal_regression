function plotKernels(temporalKernels,sessIDs,timeBack,timeForward,eventNames,params,saveLoc,fname,shuffKernels,CIFlag)

if nargin < 10
    CIFlag = 0;
end

%% Plot all events on the same axis by session

f=figure('Units','inches','Position',[4.5000 0.4236  6.031 4.4551]);
t=tiledlayout(numel(sessIDs),2);
t.TileSpacing = 'compact';
plotMin = [];
plotMax = [];
for nr = 1:numel(params.regions)
    counter = 1; 
    for ns = sessIDs
        p(nr,counter) = nexttile(t,tilenum(t,counter,nr)); hold on
        eval(sprintf('escapeMed = median(params.escapeLatency.%s{ns});',params.regions{nr})) % find median cross times
        eval(sprintf('avoidMed = median(params.avoidLatency.%s{ns});',params.regions{nr}))
        for ne = 1:numel(eventNames)
            % Plotting parameters based on which event
            if contains(eventNames{ne},'Escape')
                cmap = params.escapeC;
            elseif contains(eventNames{ne},'Shock')
                cmap = params.shockC;
            else
                cmap = params.avoidC;
            end
            if contains(eventNames{ne},'Cross')
                lstyle = '--';
            else
                lstyle = '-';
            end

            % Extract average kernel and error 
            thisKernel = eval(sprintf('temporalKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
            mu  = mean(thisKernel,1,'omitnan');
            if params.numShuff>0
                % if bootstrapping, error as standard deviation or 95% confidence interval of the iterations
                thisKernel = eval(sprintf('shuffKernels.%s.%s(:,:,ns);',params.regions{nr},eventNames{ne}));
                if CIFlag 
                    errNeg = prctile(thisKernel,2.5,1);
                    errPos = prctile(thisKernel,97.5,1);
                else
                    err = std(thisKernel,[],1);
                end
            else
                err = nansem(thisKernel,1);
            end
            if contains(eventNames{ne},'Cue')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2));
            elseif contains(eventNames{ne}, 'Shock')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + 5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Escape')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + escapeMed +5;
            elseif contains(eventNames{ne},'Cross') && contains(eventNames{ne}, 'Avoid')
                x = linspace(-timeBack(ne),timeForward(ne),size(thisKernel,2)) + avoidMed;
            end

            plot(x,mu,'Color',cmap,'LineWidth',1.5,'LineStyle',lstyle);
            if CIFlag
                patch([x flip(x)],[errNeg flip(errPos)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            else
                patch([x flip(x)],[mu-err flip(mu+err)],cmap,'FaceAlpha',.1,'EdgeColor','none')
            end
        end
        axis tight
        plotMin = cat(1,plotMin,p(nr,counter).YLim(1));
        plotMax = cat(1,plotMax,p(nr,counter).YLim(2));
        if ns == 1
            title(params.regions{nr})
        end
        if ns == max(sessIDs)
            xlabel('Time from cue on (s)')
        end
        %ylabel(sprintf('Day %s', num2str(ns)))
        counter = counter+1; % plot counter
    end
end
ylabel(t,'Kernel value')

plotMin = min(plotMin);
plotMax = max(plotMax);
for nr = 1:size(p,1)
    for ns = 1:size(p,2)
        eval(sprintf('escapeMed = median(params.escapeLatency.%s{ns});',params.regions{nr})) % find median cross times 
        eval(sprintf('avoidMed = median(params.avoidLatency.%s{ns});',params.regions{nr}))
        set(p(nr,ns),'YLim', [plotMin plotMax])
        plot(p(nr,ns),[0 0], [plotMin plotMax],'--k');
        l(3)=plot(p(nr,ns),[5 5],[plotMin plotMax],'--','Color',params.shockC);
        l(2)=plot(p(nr,ns),[escapeMed+5 escapeMed+5], [plotMin plotMax],'--','Color',[params.escapeC .5]);
        l(1)=plot(p(nr,ns),[avoidMed avoidMed], [plotMin plotMax],'--','Color',[params.avoidC .5]);
    end
end
%legend(p(2,end),l,{'Avoid cross'; 'Escape cross';'Shock'},'Box','off','Location','southeast');

exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)))




