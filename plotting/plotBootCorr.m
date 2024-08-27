function stats = plotBootCorr(realStats,shuffStats,params,eventNames,CIFlag,saveLoc,fname)

al = .05./(numel(params.sessIDs)*numel(params.regions));

%% Plot correlation coefficient for the full model 

for nr = 1:numel(params.regions)
    for ns = 1:numel(params.sessIDs)
        mu(nr,ns) = realStats{nr}{ns}.corr.full;
        clear thisShuff
        for nss = 1:params.numShuff
            thisShuff(nss) = shuffStats{nr}{ns}{nss}.corr.full;
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
p(1)=scatter(x,mu(1,:),10,'MarkerFaceColor','k','MarkerEdgeColor','none');
if CIFlag
    plot([x; x], [errNeg(1,:);errPos(1,:)],'k')
    plot([x; x], [errNeg(2,:);errPos(2,:)],'k')
else
    errorbar(x,mu(1,:),errNeg(1,:),errPos(1,:),'LineStyle','none','CapSize',0,'Color','k')
    errorbar(x,mu(2,:),errNeg(2,:),errPos(2,:),'LineStyle','none','CapSize',0,'Color','k')
end
p(2)=scatter(x,mu(2,:),10,'MarkerEdgeColor','none','MarkerFaceColor',[34 194 227]./255);
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
al = .05./(numel(params.sessIDs)*numel(params.regions) + numel(params.regions)*6);
for ne = 1:numel(eventNames)
    for nr = 1:numel(params.regions)
        for ns = 1:numel(params.sessIDs)
            %mu(nr,ns)  = realStats{nr}{ns}.corr.full-eval(sprintf('realStats{nr}{ns}.corr.%s',eventNames{ne}));
            clear thisShuff
            for nss = 1:params.numShuff
                thisShuff_full(nss) = shuffStats{nr}{ns}{nss}.corr.full;
                thisShuff_reduced(nss)  =  eval(sprintf('shuffStats{nr}{ns}{nss}.corr.%s',eventNames{ne}));
            end
            thisShuff = thisShuff_full-thisShuff_reduced;
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


    % Compare Day 1 to the rest of the days
    for nr = 1:numel(params.regions)
        for ns = 2:numel(params.sessIDs)
            thisShuff = shuff{nr,1}-shuff{nr,ns};
            CI = prctile(thisShuff,[100*al/2,100*(1-al/2)]);
            H_day(nr,ns-1) =  CI(1)>0 | CI(2)<0;
            P_day(nr,ns-1) = mean(thisShuff<=0)*2;
            %figure(); histogram(thisShuff);
        end
    end

    f=figure('Units','inches','Position',[4.5729 2.8854 2.1357 1.6253]); hold on
    x = 1:numel(params.sessIDs);
    if CIFlag
        plot([x+.1; x+.1], [errNeg(1,:);errPos(1,:)],'k')
        plot([x-.1; x-.1], [errNeg(2,:);errPos(2,:)],'k')
        y(1,:) = errPos(1,:)+.002;
        y(2,:) = errPos(2,:)+.002;
    else
        errorbar(x+.1,mu(1,:),errNeg(1,:),errPos(1,:),'LineStyle','none','CapSize',0,'Color','k')
        errorbar(x-.1,mu(2,:),errNeg(2,:),errPos(2,:),'LineStyle','none','CapSize',0,'Color','k')
        y(1,:) = mu(1,:)+errPos(1,:)+.002;
        y(2,:) = mu(2,:)+errPos(2,:)+.002;
    end
    p(1)=scatter(x+.1,mu(1,:),10,'MarkerFaceColor','k','MarkerEdgeColor','none');
    text(x(H(1,:))+.1,y(1,H(1,:)),'*','FontSize',10)
    p(2)=scatter(x-.1,mu(2,:),10,'MarkerEdgeColor','none','MarkerFaceColor',[34 194 227]./255);
    text(x(H(2,:))-.1,y(2,H(2,:)),'*','FontSize',10)
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

