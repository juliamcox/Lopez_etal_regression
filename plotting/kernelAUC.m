function kernelAUC(temporalKernels,shuffKernels,params,whichTests,whichRegions,sessIDs,CIFlag,saveLoc,fname)



[r,c] = find(whichTests); % find which events to compare
al = .05/(numel(r)*numel(sessIDs));

for nr = 1:numel(whichRegions)

    for ne = 1:numel(params.eventNames)
        for ns = sessIDs
            auc(ne,ns) = trapz(eval(sprintf('temporalKernels.%s.%s(:,1+params.timeBack(ne)*params.newFs:end,ns)''',whichRegions{nr},params.eventNames{ne})));
            % Area under the curve from bootstrapping
            for nss = 1:params.numShuff
                auc_shuff{ne}(ns,nss) = trapz(eval(sprintf('shuffKernels.%s.%s(nss,1+params.timeBack(ne)*params.newFs:end,ns)''',whichRegions{nr},params.eventNames{ne})));
            end
            auc(ne,ns)  =median(auc_shuff{ne}(ns,:));
            % Calculate error from bootstrapping
            if CIFlag
                errNeg(ne,ns) = prctile(auc_shuff{ne}(ns,:),2.5);
                errPos(ne,ns) = prctile(auc_shuff{ne}(ns,:),97.5);
            else
                errNeg(ne,ns) = std(auc_shuff{ne}(ns,:));
                errPos(ne,ns) = std(auc_shuff{ne}(ns,:));
            end
        end
    end

    x = sessIDs;
    for nrr = 1:numel(r)
        % Determine whether AUC for the events are significantly different
        for ns = sessIDs
            aucDiff  = abs(auc_shuff{r(nrr)}(ns,:))-abs(auc_shuff{c(nrr)}(ns,:));
            CI = prctile(aucDiff,[100*al/2,100*(1-al/2)]);
            H(nrr,ns) = CI(1)>0 | CI(2)<0;
        end
        f = figure('Units','inches'); f.Position(3) = 1.9194; f.Position(4) = 1.6253; hold on
        p(1) = scatter(x-.1,auc(r(nrr),:),10,'MarkerFaceColor',params.avoidC,'MarkerEdgeColor','none');
        p(2) = scatter(x+.1,auc(c(nrr),:),10,'MarkerFaceColor',params.escapeC,'MarkerEdgeColor','none');
        if CIFlag
            plot([x-.1;x-.1],[errNeg(r(nrr),:); errPos(r(nrr),:)],'Color',p(1).MarkerFaceColor)
            plot([x+.1;x+.1],[errNeg(c(nrr),:); errPos(c(nrr),:)],'Color',p(2).MarkerFaceColor)
        else
            errorbar(x-.1, auc(r(nrr),:),errNeg(r(nrr),:),errPos(r(nrr),:), 'Color',p(1).MarkerFaceColor,'LineStyle','none','CapSize',0,'LineWidth',1)
            errorbar(x+.1, auc(c(nrr),:),errNeg(c(nrr),:),errPos(c(nrr),:), 'Color',p(2).MarkerFaceColor,'LineStyle','none','CapSize',0,'LineWidth',1)
        end
        g=gca;
        text(x(H(nrr,:)),ones(sum(H(nrr,:)),1).*g.YLim(2)-.02,'*','FontSize',11)
        set(gca,'XLim',[.5 x(end)+.5])
        %legend(p,{params.eventNames{r(nrr)};params.eventNames{c(nrr)}},'Box', 'off')
        ylabel('Area under the curve')
        xlabel('Day')
        box off
        title(whichRegions{nr})
        exportgraphics(f,fullfile(saveLoc,sprintf('%s.pdf',fname)),'Append',true)

    end
end

   