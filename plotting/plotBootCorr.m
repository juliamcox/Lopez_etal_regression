function stats = plotBootCorr(realStats,shuffStats,params,eventNames,CIFlag,saveLoc,fname)

al = .05./(numel(params.sessIDs)*numel(params.regions));
plotColors = cat(1,[34 194 227]./255,[0 0 0]);

%% Plot correlation coefficient for the full model 

for nr = 1:numel(params.regions)
    for ns = params.sessIDs
        mu(nr,ns) = realStats{nr}{ns}.corr.full;
        clear thisShuff
        try
        for nss = 1:params.numShuff
            thisShuff(nss) = shuffStats{nr}{ns}{nss}.corr.full;
        end
        catch
            thisShuff = nan(1,params.numShuff);
        end
       mu(nr,ns) = median(thisShuff, 'omitmissing');
        if CIFlag
            errNeg(nr,ns) = prctile(thisShuff,2.5);
            errPos(nr,ns) = prctile(thisShuff,97.5);
        else
            errNeg(nr,ns) = std(thisShuff);
            errPos(nr,ns) = std(thisShuff);
        end
    end
end

% Fit a line to each region
x = 1:numel(params.sessIDs);
for nr = 1:numel(params.regions)
    mdl = fitlm(x',mu(nr,:)');
    y_hat(nr,:) = mdl.Coefficients.Estimate(2).*x + mdl.Coefficients.Estimate(1);
    stats.line.slope(nr)   = mdl.Coefficients.Estimate(2);
    stats.line.pval(nr,:)  = mdl.Coefficients.pValue(2);
end


f=figure('Units','inches','Position',[4.5729 2.8854 1.9427 1.6253]); hold on
plot(x,y_hat,'k')
for nr = 1:numel(params.regions)

    if CIFlag
        plot([x; x], [errNeg(nr,:);errPos(nr,:)],'k')

    else
        errorbar(x,mu(nr,:),errNeg(nr,:),errPos(nr,:),'LineStyle','none','CapSize',0,'Color',plotColors(nr,:))
    end

    p(nr)=scatter(x,mu(nr,:),10,'MarkerFaceColor',plotColors(nr,:),'MarkerEdgeColor','none');

end
g=gca;
g.YLim(1)=0;
g.XLim = [.5 x(end)+.5];
g.XTick = params.sessIDs;
g.FontSize = 7.5623;
ylabel('Correlation coefficient','FontSize',7.5623)
xlabel('Day','FontSize',7.5623)
title('Full model','FontSize',8.6426)
legend(p,params.regions,'Box','off','Location','southwest')

exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

%% Plot correlation coefficient for the full model - reduced model
al = .05./(numel(params.sessIDs)*numel(params.regions) + numel(params.regions)*numel(eventNames));
clear mu H P 
for ne = 1:numel(eventNames)
    for nr = 1:numel(params.regions)
        for ns = params.sessIDs
            %mu(nr,ns)  = realStats{nr}{ns}.corr.full-eval(sprintf('realStats{nr}{ns}.corr.%s',eventNames{ne}));
            clear thisShuff
            try
                for nss = 1:params.numShuff
                    thisShuff_full(nss) = shuffStats{nr}{ns}{nss}.corr.full;
                    thisShuff_reduced(nss)  =  eval(sprintf('shuffStats{nr}{ns}{nss}.corr.%s',eventNames{ne}));
                end
                thisShuff = thisShuff_full-thisShuff_reduced;
            catch
                thisShuff = nan(1,params.numShuff);
            end
            mu(nr,ns) = median(thisShuff,'omitmissing');
            shuff{nr,ns} = thisShuff;
            CI = prctile(thisShuff,[100*al/2,100*(1-al/2)]);
            H(nr,ns) = CI(1)>0 | CI(2)<0;
            P(nr,ns) = mean(thisShuff<=0)*2;
            if CIFlag
                errNeg(nr,ns) = prctile(thisShuff,2.5);
                errPos(nr,ns) = prctile(thisShuff,97.5);
            else
                errNeg(nr,ns) = std(thisShuff);
                errPos(nr,ns) = std(thisShuff);
            end
        end
    end

    eval(sprintf('stats.%s = P;',eventNames{ne}));


    % % Compare Day 1 to the rest of the days
    % for nr = 1:numel(params.regions)
    %     for ns = 2:numel(params.sessIDs)
    %         thisShuff = shuff{nr,1}-shuff{nr,ns};
    %         CI = prctile(thisShuff,[100*al/2,100*(1-al/2)]);
    %         H_day(nr,ns-1) =  CI(1)>0 | CI(2)<0;
    %         P_day(nr,ns-1) = mean(thisShuff<=0)*2;
    %         %figure(); histogram(thisShuff);
    %     end
    % end

    f=figure('Units','inches','Position',[4.5729 2.8854 2.1357 1.6253]); hold on
    x = (1:numel(params.sessIDs))-.2;
    for nr = 1:numel(params.regions)
        x = x+.4*(nr-1);
        if CIFlag
            plot([x; x], [errNeg(nr,:);errPos(nr,:)],'Color',plotColors(nr,:))
            y(nr,:) = errPos(nr,:)+.002;
        else
            errorbar(x,mu(1,:),errNeg(nr,:),errPos(nr,:),'LineStyle','none','CapSize',0,'Color',plotColors(nr,:))
            y(nr,:) = mu(nr,:)+errPos(nr,:)+.002;
        end
        p(nr)=scatter(x,mu(nr,:),10,'MarkerFaceColor',plotColors(nr,:),'MarkerEdgeColor','none');
        text(x(H(nr,:)),y(nr,H(nr,:)),'*','FontSize',10)
    end
    plot([.5 x(end)+.5],[0 0],':','Color',[0 0 0 .5])
    g=gca;
    g.XLim = [.5 x(end)+.5];
    g.XTick = params.sessIDs;
    g.FontSize = 7.5623;
    ylabel(sprintf('Correlation coefficient\n(full - reduced)'))
    xlabel('Day')
    title(eventNames{ne}, 'FontSize',8.6426)
    legend(p,params.regions,'Box','off')
    exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

end

%% Cacluate p value for speed predictor
if params.speedFlag
    speedEvents = fieldnames(shuffStats{1}{1}{1}.corr);
    speedEvents = speedEvents(contains(speedEvents,'speed'));

    for ne = 1:numel(speedEvents)
        x = (1:numel(params.sessIDs))-.2;
        for nr = 1:numel(params.regions)
            x = x+.4*(nr-1);
            for ns = params.sessIDs
                clear thisShuff
                try
                    for nss = 1:params.numShuff
                        thisShuff_full(nss) = shuffStats{nr}{ns}{nss}.corr.full;
                        thisShuff_reduced(nss)  =  eval(sprintf('shuffStats{nr}{ns}{nss}.corr.%s',speedEvents{ne}));
                    end
                    thisShuff = thisShuff_full-thisShuff_reduced;
                catch
                    thisShuff = nan(1,params.numShuff);
                end
                mu(nr,ns) = median(thisShuff,'omitmissing');
                shuff{nr,ns} = thisShuff;
                CI = prctile(thisShuff,[100*al/2,100*(1-al/2)]);
                H(nr,ns) = CI(1)>0 | CI(2)<0;
                P(nr,ns) = mean(thisShuff<=0)*2;
                if CIFlag
                    errNeg(nr,ns) = prctile(thisShuff,2.5);
                    errPos(nr,ns) = prctile(thisShuff,97.5);
                else
                    errNeg(nr,ns) = std(thisShuff);
                    errPos(nr,ns) = std(thisShuff);
                end
            end
        end
        
        eval(sprintf('stats.speed%s = P;',num2str(ne)));

        f=figure('Units','inches','Position',[4.5729 2.8854 2.1357 1.6253]); hold on
        x = 1:numel(params.sessIDs);
        for nr = 1:numel(params.regions)
            if CIFlag
                plot([x; x], [errNeg(nr,:);errPos(nr,:)],'Color',plotColors(nr,:))
                y(nr,:) = errPos(nr,:)+.002;
            else
                errorbar(x,mu(1,:),errNeg(nr,:),errPos(nr,:),'LineStyle','none','CapSize',0,'Color',plotColors(nr,:))
                y(nr,:) = mu(nr,:)+errPos(nr,:)+.002;
            end
            p(nr)=scatter(x,mu(nr,:),10,'MarkerFaceColor',plotColors(nr,:),'MarkerEdgeColor','none');
            text(x(H(nr,:)),y(nr,H(nr,:)),'*','FontSize',10)
        end
        plot([.5 x(end)+.5],[0 0],':','Color',[0 0 0 .5])
        g=gca;
        g.XLim = [.5 x(end)+.5];
        g.XTick = params.sessIDs;
        g.FontSize = 7.5623;
        ylabel(sprintf('Correlation coefficient\n(full - reduced)'))
        xlabel('Day')
        title(speedEvents{ne}(1:end-1), 'FontSize',8.6426)
        legend(p,params.regions,'Box','off')
        exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)
    end
end

