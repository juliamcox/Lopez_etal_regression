function stats = multCompareCorr(realStats,params,saveLoc,fname)


% Extract correlation coefficients 
for nr = 1:numel(params.regions)
    for na = 1:size(realStats{nr},1)
        for ns = 1:size(realStats{nr},2)
            corr_full(na,ns,nr) = realStats{nr}{na,ns}.corr.full;
            z_full(na,ns,nr)    = atanh(corr_full(na,ns,nr));
            if isnan(corr_full(na,ns,nr))
                for ne = 1:numel(params.eventNames)
                    corr_event{ne}(na,ns,nr) = NaN;
                    z_event{ne}(na,ns,nr) = NaN;
                end
            else
                for ne = 1:numel(params.eventNames)
                    corr_event{ne}(na,ns,nr) = eval(sprintf('realStats{nr}{na,ns}.corr.%s',params.eventNames{ne}));
                    z_event{ne}(na,ns,nr) = atanh(corr_event{ne}(na,ns,nr));
                end
            end
        end
    end
end


figure('Position',[426 63 362 686]); hold on
t = tiledlayout(numel(params.eventNames)+1,numel(params.regions));
t.Padding = "compact";
for nr = 1:numel(params.regions)
    thisCorr = corr_full(1:size(realStats{nr},1),:,nr);
    thisZ    = array2table(z_full(1:size(realStats{nr},1),:,nr));
    Time     = table([params.sessIDs]','VariableNames',{'Time'});
    stats.rm_full{nr} = ranova(fitrm(thisZ,'Var1-Var7~1',WithinDesign=Time));
    
    p(nr,1) = nexttile; hold on
    if nr == 1
        cmap = 'k';
    else
        cmap = [34 194 227]./255;
    end
    scatter(params.sessIDs,mean(thisCorr,1,'omitnan'),20,'MarkerFaceColor',cmap,'MarkerEdgeColor','none');
    errorbar(params.sessIDs,mean(thisCorr,1,'omitnan'),nansem(thisCorr,1),'k','LineStyle','none','CapSize',0)
    if stats.rm_full{nr}.pValue(1)<0.05
        plot([1 params.sessIDs(end)], [.6 .6], 'k')
        text(round(numel(params.sessIDs)/2),.65,'*')
    end
    title('Full model')
    set(gca,'XLim',[0 params.sessIDs(end)+1],'XTick',params.sessIDs,'YLim',[0 .7])
    if nr == 1
        ylabel('Correlation coefficient','FontSize',10)
    end
end

plotLims = [.08,.08,.08,.2,.1];
for ne = 1:numel(params.eventNames)
    for nr = 1:numel(params.regions)
    p(nr,1+ne) = nexttile; hold on
    thisCorr = corr_full(1:size(realStats{nr},1),:,nr)-corr_event{ne}(1:size(realStats{nr},1),:,nr);
    thisZ    = array2table(z_full(1:size(realStats{nr},1),:,nr)-z_event{ne}(1:size(realStats{nr},1),:,nr));
    Time     = table([params.sessIDs]','VariableNames',{'Time'});
    stats.rm_event{nr,ne} = ranova(fitrm(thisZ,'Var1-Var7~1',WithinDesign=Time));
    if nr == 1
        cmap = 'k';
    else
        cmap = [34 194 227]./255;
    end
    scatter(params.sessIDs,mean(thisCorr,1,'omitnan'),20,'MarkerFaceColor',cmap,'MarkerEdgeColor','none');
    errorbar(params.sessIDs,mean(thisCorr,1,'omitnan'),nansem(thisCorr,1),'k','LineStyle','none','CapSize',0)
    
    if stats.rm_event{nr,ne}.pValue(1)<0.05
        plot([1 params.sessIDs(end)], [plotLims(ne)-.02 plotLims(ne)-.02], 'k')
        text(round(numel(params.sessIDs)/2),plotLims(ne)-.01,'*')
    end
    title(params.eventNames{ne})
    set(gca,'XLim',[0 params.sessIDs(end)+1],'XTick',params.sessIDs,'YLim',[0 plotLims(ne)])
    end
end
xlabel(t,'Day')
ylabel(t,'Correlation coefficient (full-reduced)       ')



