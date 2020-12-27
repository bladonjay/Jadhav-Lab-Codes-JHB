% SummaryObjPlaceInteractions
%{
 this is the summary figures for the interaction between object coding
 and place coding.  The central question is whether object coding or
 responsive cells show any difference in their spatial tendencies.

basically the research and questions are as follows:

---------------------------
Taxidis
 the taxidis result is that across epochs place cells are more stable than odor cells


-------------------------------------
Terada
the terada result is that odor cells will phase precess just like place cells


--------------------------
the li xyu liu liu zhu zhang etcx. a distinct entorhinal circut paper
main result was that the spatial properties of odor cells are identical
to non odor cells
-this paper really didnt focus on this result


------------------------------
Igarashi:

they dont address the spatial characteristics, but...

-PV coding after cue sampling should be stronger when they know the task better
-unit selectivity, coherence, and performance are all correlated


-----------------------------
Fujisawa paper:
-PFC and VTA neurons cohere to 4-hz oscillations
-outcome predicting pfc cells cohere to 4 hz oscillation
-i think these are pfc splitter cells

-----------------------------
Rangel paper:
-theta, beta, and gamma coherence should be better in INS and pyrs for
    correct trials, its not related to firing rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reasonable questions:

1. do the odor fields look anything like place fields?
    -cat single trial peaks, how wide are the fields?
    -do odor fields theta precess even though theta is weak? (maybe rr
    precess)
2. are odor cells equally likely to have place fields on the maze? and are
    place cells more likely to be odor responsive?
3. do cells that respond to odors have different shaped place fields?
4. do odor selectivce cells have same or different route selectivity
5. are odor cells more stable or less stable than odor selective cels
6. are splitter cells more likely to be modulated by odor?

Answered questions:

-do object responsive or object selective cells have a different propensity
to have place fields?
-do these cells have place fields with a different peak, width, or
information content?
-Are these cells more apt to select for one trajectory over another in the
same place?  e.g. are they splitters or are they directionally sensitive?

%}


%%
% get a superrat of units
SuperUnits=orderfields(SuperRat(1).units);

for i=2:length(SuperRat)
    theseunits=orderfields(SuperRat(i).units);
    % for some reason csi is in some of these structs
    if isfield(theseunits,'csi'), theseunits=rmfield(theseunits,'csi'); end
    SuperUnits=[SuperUnits theseunits];
end
fprintf('Number of well isolated units %d \n',...
    sum(cellfun(@(a) ~contains(a,'mua','IgnoreCase',true), {SuperUnits.type})))

hcfilt=cellfun(@(a) contains(a,'CA1'), {SuperUnits.area});
fprintf('CA1 %d pyrs, %d IN''s \n', sum(cellfun(@(a) contains(a,'pyr'), {SuperUnits(hcfilt).type})),...
    sum(cellfun(@(a) contains(a,'in'), {SuperUnits(hcfilt).type})));

pfcfilt=cellfun(@(a) contains(a,'PFC'), {SuperUnits.area});
fprintf('pfc %d pyrs, %d IN''s \n', sum(cellfun(@(a) contains(a,'pyr'), {SuperUnits(pfcfilt).type})),...
    sum(cellfun(@(a) contains(a,'in'), {SuperUnits(pfcfilt).type})));
%%
% 1. do cells that respond to odors respond to place differently
% First, how many of these cells spike during the outbound run bouts
region={'CA1','PFC'};
type='pyr';
for r=1:2
    % the only pool in which it does matter whether
    cellfilt=cellfun(@(a) contains(a,region{r},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type});
    
    cellPool=SuperUnits(cellfilt);
    
    % what %% of cells overall have a pf
    
    % 1 if left preferring, -1 if right preferring, 0 if no pref
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) a(1)>0, {cellPool.OdorSelective})-.5)*2);
    % -1 if downreg, 0 if no reg, 1 if upreg (This ough tto work as only one
    % cell is going to have opposite firing on the runs)
    objRes=cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}).*...
        ((cellfun(@(a) any(a((a(:,4)==min(a(:,4))),2)>0),{cellPool.OdorResponsive})-.5)*2);
    
    
    % ok %% of cells with 1 and %% of cells with more than 1 pf that arent obj
    % for forward and reverse or just foreward?
    PFcounts=cell2mat(cellfun(@(a) a(1:2), {cellPool.FiresDuringRun},'UniformOutput',false)');
    
    bary=[nanmean(sum(PFcounts(objSel==0 & objRes==0,:),2)==1)...
        nanmean(sum(PFcounts(objSel==0 & objRes~=0,:),2)==1)...
        nanmean(sum(PFcounts(objSel~=0,:),2)==1);...
        nanmean(sum(PFcounts(objSel==0 & objRes==0,:),2)>1)...
        nanmean(sum(PFcounts(objSel==0 & objRes~=0,:),2)>1)...
        nanmean(sum(PFcounts(objSel~=0,:),2)>1)]*100;
    figure;
    b=bar(bary','stacked');
    set(gca,'XTickLabel',{sprintf('Non Odor Responsive (%d)',sum(~objSel & objRes==0)),...
        sprintf('Odor Responsive (%d)', sum(objRes~=0 & ~objSel)),...
        sprintf('Odor Selective (%d)', sum(objSel))});
    legend('Spikes on one rounte','on multiple routes');
    ylabel('% of units with spikes');
    xtips = b(2).XEndPoints; ytips = b(1).YEndPoints;  ytips2 = b(2).YEndPoints;
    labels = {sum(sum(PFcounts(objSel==0 & objRes==0,:),2)==1) sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)==1) sum(sum(PFcounts(objSel~=0,:),2)==1)};
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    labels2 = {sum(sum(PFcounts(objSel==0 & objRes==0,:),2)>1) sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)>1) sum(sum(PFcounts(objSel~=0,:),2)>1)};
    text(xtips,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    ylim([0 100]);
    sgtitle(sprintf('Outbound Run Firing of %s cells from %s, out of a pool of %d', type, region{r}, length(cellPool)));
    xtickangle(340); box off;
    
    fprintf('%s numbers: Nonresponsive tot %d: one %d (%.f): two %d (%.f) \n',...
        region{r},sum(objSel==0 & objRes==0), labels{1}, bary(1,1), labels2{1},bary(2));
    fprintf('%s numbers: Responsive tot %d: one %d (%.f): two %d (%.f) \n',...
        region{r},sum(objSel==0 & objRes~=0), labels{2}, bary(1,2), labels2{1},bary(2,2));
    fprintf('%s numbers: Selective tot %d: one %d (%.f): two %d (%.f) \n',...
        region{r},sum(objSel~=0), labels{3}, bary(1,3), labels2{3},bary(1,3));
    
    
    propsbig=[sum(sum(PFcounts(objSel==0 & objRes==0,:),2)==0)...
        sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)==0)...
        sum(sum(PFcounts(objSel~=0,:),2)==0);...
        sum(sum(PFcounts(objSel==0 & objRes==0,:),2)==1)...
        sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)==1)...
        sum(sum(PFcounts(objSel~=0,:),2)==1);...
        sum(sum(PFcounts(objSel==0 & objRes==0,:),2)>1)...
        sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)>1)...
        sum(sum(PFcounts(objSel~=0,:),2)>1)];
    [a,b,c]=chi2indep(propsbig);
    fprintf('%s Chi2 test, # routs with spikes across odor repsonses chi2(%d)= %.2f p=%.2e \n',...
        region{r},c,b,a)
end
%% okay now to plot odor period rates vs run rates
 % no filters on minimum rates

odormeans=cellfun(@(a) nanmean(a), {SuperUnits.OdorMeans});
runmeans=cellfun(@(a) mean(cell2mat(a)), {SuperUnits.RunRates});
figure;
clear sp si; mycolors=lines(2);
for i=1:length(region)
    regionfilt=cellfun(@(a) contains(a,region{i}), {SuperUnits.area});
    PYRfilt=cellfun(@(a) contains(a,'pyr'), {SuperUnits.type});
    INfilt=cellfun(@(a) contains(a,'in'), {SuperUnits.type});
    nanfilt=~isnan(odormeans) & ~isnan(runmeans);
    
    [a,b]=corr(odormeans(regionfilt & PYRfilt & nanfilt)',runmeans(regionfilt & PYRfilt & nanfilt)','type','spearman');
    [a2,b2]=corr(odormeans(regionfilt & INfilt & nanfilt)',runmeans(regionfilt & INfilt & nanfilt)','type','spearman');
    mdl = fitlm(odormeans(regionfilt & PYRfilt & nanfilt)',runmeans(regionfilt & PYRfilt & nanfilt)');
    
    fprintf('%s pyr odor rate to run rate n=%d corr=%.2f p=%.2e, %s in corr=%.2f p=%.2e\n',...
        region{i},sum(regionfilt & PYRfilt & nanfilt),a,b,region{i},a2,b2)
    
    sp(i)=subplot(2,2,i);
    %xvals=log(odormeans(regionfilt & PYRfilt & nanfilt));
    %yvals=log(runmeans(regionfilt & PYRfilt & nanfilt));
    
    xvals=odormeans(regionfilt & PYRfilt & nanfilt);
    yvals=runmeans(regionfilt & PYRfilt & nanfilt);
    
    okuse=~isinf(abs(xvals)) & ~isinf(abs(yvals)) & yvals~=0 & xvals~=0;
    scatter(xvals,yvals,4,mycolors(1,:),'filled');
    xlabel('Rate during odorpoke'); ylabel('Rate during run');  
    title(sprintf('%s r2=%.2f p=%.2e',region{i},a,b)); legend('pyr','fit');
    hold on;
    [r,m,b]=regression(log(xvals(okuse)),log(yvals(okuse)));
    
    set(gca,'Xscale','log','Yscale','log');
    %set(gca,'XTick',[-4:2:4],'XTickLabel',{'10^-^4','10^-^2','1','10^2','10^4'});
   % set(gca,'YTick',[-6:2:4],'YTickLabel',{'10^-^6','10^-^4','10^-^2','1','10^2','10^4'});
    xlims=log(get(gca,'XLim'));
    plot(exp([xlims])',exp([[xlims]*m+b]'),'r'); set(gca,'Xlim',xlims);
    set(gca,'XLim',[0.1 100]);
    
    si(i)=subplot(2,2,i+2);
    scatter(odormeans(regionfilt & INfilt & nanfilt),runmeans(regionfilt & INfilt & nanfilt),4,mycolors(2,:),'filled');
    xlabel('Rate during odorpoke'); ylabel('Rate during run');set(gca,'Xscale','log','Yscale','log');
    title(sprintf('%s r2=%.2f p=%.2e',region{i},a2,b2)); legend('IN');
    
    
end
%linkaxes(sp);
%linkaxes(si);


%% does place field prevalence differ across odor responses (SELECTIVE AND RESPONSIVE)
% now filter out the cells who dont spike during runs
fprintf('\n \n');
% I think these cells will be called 'task relevant firing' e.g. fired on
% tyhe run, or fired during the odor sampling epoch (at least N spikes)
for i=1:length(region)
    
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) any(a),{SuperUnits.FiresDuringRun});
   
    
    cellPool=SuperUnits(cellfilt);
    % what %% of cells that fire on a run have a pf on that run
    % 1 if left preferring, -1 if right preferring, 0 if no pref
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) a(1)>0, {cellPool.OdorSelective})-.5)*2);
    % -1 if downreg, 0 if no reg, 1 if upreg
    % to get direction, use the min p index (they're always the same direction)
    
    objRes=cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}).*...
        ((cellfun(@(a) any(a((a(:,4)==min(a(:,4))),2)>0),{cellPool.OdorResponsive})-.5)*2); % we cant do
    
    
    % condense into cell classes
    cellclasses=double(objSel==0 & objRes==0);
    cellclasses(objSel==0 & objRes~=0)=2;
    cellclasses(objSel~=0)=3;
    
    % ok %% of cells with 1 and %% of cells with more than 1 pf that arent obj
    PFcounts=cell2mat(cellfun(@(a) a(1:2), {cellPool.PFexist},'UniformOutput',false)');
    
    nbars=3;
    if nbars==3
   bary=[nanmean(sum(PFcounts(cellclasses==1,:),2)==1)...          % nonselective
        nanmean(sum(PFcounts(cellclasses==2,:),2)==1)...            % odor responsive
        nanmean(sum(PFcounts(cellclasses==3,:),2)==1);...           % odor selective
        nanmean(sum(PFcounts(cellclasses==1,:),2)>1)...
        nanmean(sum(PFcounts(cellclasses==2,:),2)>1)...
        nanmean(sum(PFcounts(cellclasses==3,:),2)>1)]*100;
    elseif nbars==2
        
     bary=[nanmean(sum(PFcounts(cellclasses==1,:),2)==1)...          % nonselective
        nanmean(sum(PFcounts(cellclasses>1,:),2)==1);...            % odor responsive
        nanmean(sum(PFcounts(cellclasses==1,:),2)>1)...
        nanmean(sum(PFcounts(cellclasses>1,:),2)>1)]*100;
    end
    figure;
    stackbar=false;
    if stackbar
        b=bar(bary','stacked');
        legend('PF on one outbound rounte','PF on both outbounds');
        ylabel('% of units with fields');
        xtips = b(2).XEndPoints; ytips = b(1).YEndPoints;  ytips2 = b(2).YEndPoints;
        labels = {sum(sum(PFcounts(cellclasses==1,:),2)==1)...
            sum(sum(PFcounts(cellclasses==2,:),2)==1)...
            sum(sum(PFcounts(cellclasses==3,:),2)==1)};
        text(xtips,ytips,labels(1:length(xtips)),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        labels2 = {sum(sum(PFcounts(cellclasses==1,:),2)>1)...
            sum(sum(PFcounts(cellclasses==2,:),2)>1)...
        sum(sum(PFcounts(cellclasses==3,:),2)>1)};
        text(xtips,ytips2,labels2(1:length(xtips)),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
    else
        b=bar(sum(bary));
        xtips = b.XEndPoints; ytips = b.YEndPoints;
        labels1 = {sum(sum(PFcounts(cellclasses==1,:),2)>0) sum(sum(PFcounts(cellclasses==2,:),2)>0) sum(sum(PFcounts(cellclasses==3,:),2)>0)};
        text(xtips,ytips,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        legend('Spatially Selective Fragment');
    end
    mylabels={sprintf('Non Odor Responsive (%d)',sum(cellclasses==1)),...
        sprintf('Odor Responsive (%d)', sum(cellclasses==2)),...
        sprintf('Odor Selective (%d)', sum(cellclasses==3))};
    set(gca,'XTickLabel',mylabels(1:length(xtips)));
    ylim([0 100]);
    sgtitle(sprintf('Place Fields of %s cells from %s, out of a pool of %d active', type, region{i}, length(cellPool)));
    xtickangle(340); box off;
    
    
    % first row: no pfs, nonres res sel
    % second row: one pf nonres res sel
    % third row: 2pf no res     res sel
    % proportions of 1 and 2 fields
    propsbig=[sum(sum(PFcounts(cellclasses==1,:),2)==0)...        % nonselective
        sum(sum(PFcounts(cellclasses==2,:),2)==0)...            % odor responsive
        sum(sum(PFcounts(cellclasses==3,:),2)==0);...           % odor selective
        sum(sum(PFcounts(cellclasses==1,:),2)==1)...            % nonselective
        sum(sum(PFcounts(cellclasses==2,:),2)==1)...            % odor responsive
        sum(sum(PFcounts(cellclasses==3,:),2)==1);...           % odor selective
        sum(sum(PFcounts(cellclasses==1,:),2)>1)...
        sum(sum(PFcounts(cellclasses==2,:),2)>1)...
        sum(sum(PFcounts(cellclasses==3,:),2)>1)];

    % proportions of just one field
    propssmall=[sum(sum(PFcounts(cellclasses==1,:),2)==0)...        % nonselective
        sum(sum(PFcounts(cellclasses==2,:),2)==0)...            % odor responsive
        sum(sum(PFcounts(cellclasses==3,:),2)==0);...           % odor selective
        sum(sum(PFcounts(cellclasses==1,:),2)>0)...            % nonselective
        sum(sum(PFcounts(cellclasses==2,:),2)>0)...            % odor responsive
        sum(sum(PFcounts(cellclasses==3,:),2)>0)];
    [a]=1-binocdf(sum(PFcounts(cellclasses>1)),sum(cellclasses>1),...
        nanmean(PFcounts(cellclasses==1)));
    
    % the chi2 tests...
    
    fprintf('%s p(place selective) given no odor response =%.3f, (%d of %d)\n',...
        region{i},propssmall(2,1)/sum(propssmall(:,1)),propssmall(2,1),sum(propssmall(:,1)));
    
    [y]=binocdf(propssmall(2,2),sum(propssmall(:,2)),propssmall(2,1)/sum(propssmall(:,1)));
    fprintf('%s binomial test on p(spatial field |odorresponsive)=%.2f  (%d of %d cells)  pval=%.4f\n',...
        region{i},propssmall(2,2)/sum(propssmall(:,2)),propssmall(2,2),sum(propssmall(:,2)),1-y);
    
    [y]=binocdf(propssmall(2,3),sum(propssmall(:,3)),propssmall(2,1)/sum(propssmall(:,1)));
    fprintf('%s binomial test on p(spatial field |odor selective)=%.2f  (%d of %d cells)  pval=%.4f\n',...
        region{i},propssmall(2,3)/sum(propssmall(:,3)),propssmall(2,3),sum(propssmall(:,3)),1-y);
    
    fprintf('%s grand tot %d Bino Resp only, pnull=%.2f, preal=%.2f, p%.2e \n',...
        region{i},sum(cellclasses>0),nanmean(sum(PFcounts(cellclasses==1,:),2)>0),nanmean(sum(PFcounts(cellclasses>1,:),2)>0),a);
    
    
end
%%
%% does place field prevalence differ across odor responses (ONLY RESPONSIVE)
% now filter out the cells who dont spike during runs
fprintf('\n \n');
% I think these cells will be called 'task relevant firing' e.g. fired on
% tyhe run, or fired during the odor sampling epoch (at least N spikes)
for i=1:length(region)
    
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) any(a),{SuperUnits.FiresDuringRun});
   
    
    cellPool=SuperUnits(cellfilt);
    % what %% of cells that fire on a run have a pf on that run
    % 1 if left preferring, -1 if right preferring, 0 if no pref
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) a(1)>0, {cellPool.OdorSelective})-.5)*2);
    % -1 if downreg, 0 if no reg, 1 if upreg
    % to get direction, use the min p index (they're always the same direction)
    
    objRes=cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}).*...
        ((cellfun(@(a) any(a((a(:,4)==min(a(:,4))),2)>0),{cellPool.OdorResponsive})-.5)*2); % we cant do
    
    
    % condense into cell classes
    cellclasses=double(objSel==0 & objRes==0);
    cellclasses(objSel==0 & objRes~=0)=2;
    cellclasses(objSel~=0)=3;
    
    % ok %% of cells with 1 and %% of cells with more than 1 pf that arent obj
    PFcounts=cell2mat(cellfun(@(a) a(1:2), {cellPool.PFexist},'UniformOutput',false)');
    
    
    
    bary=[nanmean(sum(PFcounts(cellclasses==1,:),2)==1)...          % nonselective
        nanmean(sum(PFcounts(cellclasses>1,:),2)==1);...            % odor responsive
        nanmean(sum(PFcounts(cellclasses==1,:),2)>1)...
        nanmean(sum(PFcounts(cellclasses>1,:),2)>1)]*100;

    figure;
    stackbar=false;
    if stackbar
        b=bar(bary','stacked');
        legend('PF on one outbound rounte','PF on both outbounds');
        ylabel('% of units with fields');
        xtips = b(2).XEndPoints; ytips = b(1).YEndPoints;  ytips2 = b(2).YEndPoints;
        labels = {sum(sum(PFcounts(cellclasses==1,:),2)==1)...
            sum(sum(PFcounts(cellclasses==2,:),2)==1)...
            sum(sum(PFcounts(cellclasses==3,:),2)==1)};
        text(xtips,ytips,labels(1:length(xtips)),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        labels2 = {sum(sum(PFcounts(cellclasses==1,:),2)>1)...
            sum(sum(PFcounts(cellclasses==2,:),2)>1)...
        sum(sum(PFcounts(cellclasses==3,:),2)>1)};
        text(xtips,ytips2,labels2(1:length(xtips)),'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
    else
        b=bar(sum(bary));
        xtips = b.XEndPoints; ytips = b.YEndPoints;
        labels1 = {sum(sum(PFcounts(cellclasses==1,:),2)>0) sum(sum(PFcounts(cellclasses >1,:),2)>0)};
        text(xtips,ytips,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');
        legend('Spatially Selective Fragment');
    end
    mylabels={sprintf('Non Odor Responsive (%d)',sum(cellclasses==1)),...
        sprintf('Odor Responsive (%d)', sum(cellclasses>1))};
    set(gca,'XTickLabel',mylabels(1:length(xtips)));
    ylim([0 100]);
    sgtitle(sprintf('Place Fields of %s cells from %s, out of a pool of %d active', type, region{i}, length(cellPool)));
    xtickangle(340); box off;
    
    
    % first row: no pfs, nonres res sel
    % second row: one pf nonres res sel
    % third row: 2pf no res     res sel
    % proportions of 1 and 2 fields
   
    % proportions of just one field
    propssmall=[sum(sum(PFcounts(cellclasses==1,:),2)==0)...        % nonselective
        sum(sum(PFcounts(cellclasses>1,:),2)==0);...           % odor selective
        sum(sum(PFcounts(cellclasses==1,:),2)>0)...            % nonselective
        sum(sum(PFcounts(cellclasses>1,:),2)>0)];
    [a]=1-binocdf(sum(PFcounts(cellclasses>1)),sum(cellclasses>1),...
        nanmean(PFcounts(cellclasses==1)));
    
    % the chi2 tests...
    
    fprintf('%s p(place selective) given no odor response =%.3f, (%d of %d)\n',...
        region{i},propssmall(2,1)/sum(propssmall(:,1)),propssmall(2,1),sum(propssmall(:,1)));
    
    [y]=binocdf(propssmall(2,2),sum(propssmall(:,2)),propssmall(2,1)/sum(propssmall(:,1)));
    fprintf('%s binomial test on p(spatial field |odorresponsive)=%.2f  (%d of %d cells)  pval=%.6f\n',...
        region{i},propssmall(2,2)/sum(propssmall(:,2)),propssmall(2,2),sum(propssmall(:,2)),1-y);
    
end
%% do place field characteristics differ? e.g. mean rate, p (place selective), 
for i=1:length(region)
    
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) any(a),{SuperUnits.FiresDuringRun});
    cellPool=SuperUnits(cellfilt);
    % what %% of cells that fire on a run have a pf on that run
    % 1 if left preferring, -1 if right preferring, 0 if no pref
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) a(1)>0, {cellPool.OdorSelective})-.5)*2);
    % -1 if downreg, 0 if no reg, 1 if upreg
    % to get direction, use the min p index (they're always the same direction)
    objRes=cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}).*...
        ((cellfun(@(a) any(a((a(:,4)==min(a(:,4))),2)>0),{cellPool.OdorResponsive})-.5)*2); % we cant do
    
    % condense into cell classes (1 nonsel/resp, 2 resp, 3 sel)
    cellclasses=double(objSel==0 & objRes==0); cellclasses(objSel==0 & objRes~=0)=2; cellclasses(objSel~=0)=3;
    PFcounts=cell2mat(cellfun(@(a) a(1:2), {cellPool.PFexist},'UniformOutput',false)');
    
    % splitter (1/0)
    isSplitter=cellfun(@(a) a(1)>0 & a(2)<.05, {cellPool.SplitterScore});
    
    % this could be the overall firing rate...
    overallrate=(cellfun(@(a)  a, {cellPool.meanrate}));
    figure;
    % maybe something like this: skinny boxscatter, with all the cells with
    %sgtitle(sprintf('%s from %s',type,region{i}));
    subplot(1,3,1);
    % boxscatter all of them
    thisx=double(sum(PFcounts,2)>0)+isSplitter'*2; thisx(thisx==3)=2;
    boxScatterplot(overallrate,thisx,'XLabels',{'no pf','PF','Splitter'},...
        'Position',[],'YLabel','Mean Rate');
    xtickangle(340); box off;
    
    [a,tbl,c]=kruskalwallis(overallrate,thisx,'off');
    fprintf('%s grand mean rate spatial types: KruskalWallis Chi2(%d)=%.2f, p=%.2e \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    
    
    subplot(1,3,2);
    boxScatterplot(overallrate,cellclasses,'XLabels',{'Odor Unresponsive',...
        'Odor Responsive','Odor Selective'},'Position',[],'YLabel','Mean Rate');
    xtickangle(340); box off;
    
    
    [a,tbl,c]=kruskalwallis(overallrate,cellclasses,'off');
    fprintf('%s grand mean rate odor types: KruskalWallis Chi2(%d)=%.2f, p=%.2e \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    
    [p,~,stats]=ranksum(overallrate(cellclasses==1),overallrate(cellclasses>1));
    fprintf('%s Mean Rate: n=%d Ranksum z=%.2f, p=%.4f \n',region{i},length(overallrate),stats.zval,p);
    title(sprintf('%s Mean Rate: n=%d Ranksum z=%.2f, p=%.2e \n',region{i},length(overallrate),stats.zval,p));
    
    
    % maybe something like this: skinny boxscatter, with all the cells with
    
  
    subplot(1,3,3);
      cellclasses2=(cellclasses>1)+1;
      boxScatterplot(overallrate,cellclasses2,'XLabels',{'Odor Unresponsive',...
        'Odor Responsive'},'Position',[],'YLabel','Mean Rate');
    xtickangle(340); box off;
    
    
    [p,h,stats]=ranksum(overallrate(cellclasses2==1),overallrate(cellclasses2==2));
    fprintf('%s grand mean rate odor types: ranksum z=%.2f, p=%.2e \n',region{i},stats.zval,p);
    title(sprintf('%s Mean Rate: n=%d Ranksum z=%.2f, p=%.2e \n',...
        region{i},length(overallrate),stats.zval,p));
    
    
    
    
end
%% the inverse of this is what proportion of place cells
% have object responsiveness or object selectivity?
for i=1:length(region)
    % the pool is in this region and pyrs
    cellPool=SuperUnits(cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) any(a),{SuperUnits.FiresDuringRun}));
    
    % what %% of cells that fire on a run have a pf on that run
    
    % 1 if left preferring, -1 if right preferring, 0 if no pref
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) a(1)>0, {cellPool.OdorSelective})-.5)*2);
    % -1 if downreg, 0 if no reg, 1 if upreg
    objRes=cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}); %.*...
    %((cellfun(@(a) a(2)>0,{PFtest.OdorResponsive})-.5)*2); % we cant do
    %this right now...
    
    % ok %% of cells with 1 and %% of cells with more than 1 pf that arent obj
    % for forward and reverse or just foreward?
    PFcounts=cell2mat(cellfun(@(a) sum(a(1:2)), {cellPool.PFexist},'UniformOutput',false)');
    
    % so itll be 0 pf, 1pf, 2+pf, stacks are resp and ontop sel
    bary=[nanmean(objRes(PFcounts==0)~=0 & objSel(PFcounts==0)==0)...  % %% of 0 pf units who are selective
        nanmean(objRes(PFcounts==1)~=0 & objSel(PFcounts==1)==0)...   % %% of 1 pf units who are responsive
        nanmean(objRes(PFcounts>1)~=0 & objSel(PFcounts>1)==0);...  % %% of 2 pf units who are responsive
        nanmean(objSel(PFcounts==0)~=0)...
        nanmean(objSel(PFcounts==1)~=0)...    % %% of 1 pf units who are responsive
        nanmean(objSel(PFcounts>1)~=0)]*100;
    
    figure;
    b=bar(bary','stacked');
    set(gca,'XTickLabel',{sprintf('No Place fields (%d)',sum(PFcounts==0)),...
        sprintf('One PF (%d)', sum(PFcounts==1)),...
        sprintf('Two + fields (%d)', sum(PFcounts>1))});
    legend('Object Responsive','Object Selective');
    ylabel('% of units');
    xtips = b(2).XEndPoints; ytips = b(1).YEndPoints;  ytips2 = b(2).YEndPoints;
    labels = {sum(objRes(PFcounts==0)~=0 & objSel(PFcounts==0)==0)...
        sum(objRes(PFcounts==1)~=0 & objSel(PFcounts==1)==0)...
        sum(objRes(PFcounts>1)~=0 & objSel(PFcounts>1)==0)};
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    labels2 = {sum(objSel(PFcounts==0)~=0)...
        sum(objSel(PFcounts==1)~=0)...
        sum(objSel(PFcounts>1)~=0)};
    text(xtips,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    ylim([0 100]);
    sgtitle(sprintf('Place Fields of %s cells from %s, out of a pool of %d', type, region{i}, length(cellPool)));
    xtickangle(340); box off;
    
    propsbig=[sum(objRes(PFcounts==0)==0)...
        sum(objRes(PFcounts==1)==0)...
        sum(objRes(PFcounts>1)==0);...
        sum(objRes(PFcounts==0)~=0 & objSel(PFcounts==0)==0)...  % %% of 0 pf units who are responsive
        sum(objRes(PFcounts==1)~=0 & objSel(PFcounts==1)==0)...   % %% of 1 pf units who are responsive
        sum(objRes(PFcounts>1)~=0 & objSel(PFcounts>1)==0);...  % %% of 2 pf units who are responsive
        sum(objSel(PFcounts==0)~=0)...
        sum(objSel(PFcounts==1)~=0)...    % %% of 1 pf units who are selective
        sum(objSel(PFcounts>1)~=0)];
    
    [a,b,c]=chi2indep(propsbig);
    fprintf('%s chi2 test for obj sel given pf counts chi2(%d)=%.2f p=%e \n',region{i},c,b,a);
end
%%
% do the place fields of odor responsive/selective cells differ?
% width, in vs out of field rate, peak, and sparsity

% this only assays outbound runs, NOT INBOUND RUNS...

% each cell can have 4 pfs, but we ought to probably just take one? or just
% the outbound pfs?
% so plot histograms or cdfs or boxscatterplots
% I think we can take all the place fields of a given cell



% this plots scatterboxplots

mycolors=lines(3);
for i=1:length(region)
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type});
    
    cellPool=SuperUnits(cellfilt);
    % now I think i can probably take each unit and replicate it twice for out
    % l and out r and then cull those without place fields
    % first gather the place fields into a super long vector
    allPFs=cell2mat(cellfun(@(a) a(1:2), {cellPool.PFexist},'UniformOutput',false));
    
    % replacate this twice per cell
    % 1 if left preferring, -1 if right preferring, 0 if no pref
    objSel=repmat(cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) a(1)>0, {cellPool.OdorSelective})-.5)*2),2,1);
    %objSel=objSel';
    % -1 if downreg, 0 if no reg, 1 if upreg
    objRes=repmat(cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}),2,1); %.*...
    %((cellfun(@(a) a(2)>0,{PFtest.OdorResponsive})-.5)*2),2,1);
    %objRes=objRes';
    
    % how many cells w pfs are odor selective or responsive...
    % 0 if no pf, 1 if pf, 2 if pf and odor resp, 3 if pf odor sel
    celltypes=double(allPFs(:)==1 & objSel(:)==0 & objRes(:)==0); % 1 if it has pf, but no odor resp
    celltypes(allPFs(:)==1 & objRes(:)~=0 & objSel(:)==0)=2; % 2 if has pf and odor resp but not sel
    celltypes(allPFs(:)==1 & objSel(:)~=0)=3; % 3 if pf, and selective
    % overwrite, all 3s are now 2's for two bars (for now)
    celltypes(celltypes==3)=2;
    
    typename={'no pf','pf no odor','pf odor resp','pf odor sel'};
    typecounts=accumarray(celltypes+1,1);
    for k=1:length(typecounts)
        fprintf('Field counts for %s in %s are %d \n', typename{k}, region{i},typecounts(k));
    end
    
    figure;
    % first is pf peak rate
    PFpeaks=cell2mat(cellfun(@(a) a.PFmax(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    subplot(1,2,1);
    boxScatterplot(PFpeaks(celltypes>0),celltypes(celltypes>0),'xLabels',...
        {'Odor Unresponsive','Odor Responsive','Odor Selective'},'position',[]);
    ylabel(sprintf('%s Place Field Peak Rate (Hz)', region{i}));  xtickangle(340); box off;
    [~,tbl]=kruskalwallis(PFpeaks(celltypes>0),celltypes(celltypes>0),'off');
    fprintf('%s Peak Rate: KruskalWallis Chi2(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    [p,~,stats]=ranksum(PFpeaks(celltypes==1),PFpeaks(celltypes>1));
    fprintf('%s Peak Rate: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p);
    title(sprintf('%s Peak Rate: %d tot \n Ranksum z=%.2f, p=%.4f \n',...
        region{i},sum(celltypes>0),stats.zval,p));
    fprintf('\n');
    
    
    % this does boxscatterplots, mainly cause they make everything look the
    % Place field Width
    subplot(1,2,2);
    PFwidth=cell2mat(cellfun(@(a) a.PFsize(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    boxScatterplot(PFwidth(celltypes>0),celltypes(celltypes>0),'xLabels',...
        {'Odor Unresponsive','Odor Responsive','Odor Selective'},'Position',[]);
    ylabel('Place Field Width (% of maze)');  xtickangle(340); box off;
    
    [p,tbl,stats]=kruskalwallis(PFwidth(celltypes>0),celltypes(celltypes>0),'off');
    fprintf('%s PF Width: Kruskal Wallis Chi2(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    [p,~,stats]=ranksum(PFwidth(celltypes==1),PFwidth(celltypes>1));
    fprintf('%s PFwidth: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p);
    title(sprintf('%s PFwidth: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p));
    fprintf('\n');
    
    
    % next is sparsity
    figure;
    subplot(1,4,1);
    PFsparsity=cell2mat(cellfun(@(a) a.sparsity(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    boxScatterplot(log2(PFsparsity(celltypes>0)),celltypes(celltypes>0),'xLabels',...
        {'Odor Unresponsive','Odor Responsive','Odor Selective'},'position',[]);
    ylabel('Log_2 Place Field sparsity');  xtickangle(340); box off;
    %[p,tbl,stats]=anovan(PFsparsity(celltypes>0),celltypes(celltypes>0),'display','off');
    [p,tbl,stats]=kruskalwallis(PFsparsity(celltypes>0),celltypes(celltypes>0),'off');
    % [c,m] = multcompare(stats,'Display','off');
    fprintf('%s PF Sparsity: KruskalWallis f(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    title(sprintf('%s PF Sparsity: KruskalWallis f(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6}));
    [p,~,stats]=ranksum(PFsparsity(celltypes==1),PFsparsity(celltypes>1));
    fprintf('%s PFsparsity: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p);
    title(sprintf('%s PFsparsity: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p));
    fprintf('\n');
    
    % next is information in bits/spike
    subplot(1,4,2);
    PFinfo=cell2mat(cellfun(@(a) a.info(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    boxScatterplot(log2(PFinfo(celltypes>0)),celltypes(celltypes>0),'xLabels',...
        {'Odor Unresponsive','Odor Responsive','Odor Selective'},'position',[]);
    ylabel('Log_2 Place Field information (bits/spike)'); xtickangle(340); box off;
    
    [p,tbl,stats]=kruskalwallis(PFinfo(celltypes>0),celltypes(celltypes>0),'off');
    fprintf('%s PF Information: anovan f(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    [p,~,stats]=ranksum(PFinfo(celltypes==1),PFinfo(celltypes>1));
    fprintf('%s PFinfo: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p);
    title(sprintf('%s PFinfo: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p));
    fprintf('\n');
    
    % next is field gain (this is a hard one...)
    % here though, I am only using cells whose field isnt on the edge...
    
    % basically the pf has to be in the middle two quarters
    cumct=1;
    for j=1:length(cellPool)
        for k=1:2
            try
                linfield=cellPool(j).LinPlaceFields{1}(k,:);
                [pfmax,pfmaxpos]=max(linfield);
                if pfmaxpos>5 && pfmaxpos<95
                    pfstart=find(linfield(1:pfmaxpos) < pfmax*.25,1,'last');
                    if isempty(pfstart), pfstart=1; end
                    pfend=pfmaxpos+find(linfield(pfmaxpos:end) < pfmax*.25,1,'first');
                    if isempty(pfend), pfend=99; end
                    if ~isnan(pfend) && ~isnan(pfstart)
                        fieldgain=nanmean(linfield(pfstart+1:pfend-1))-nanmean([linfield(1:pfstart) linfield(pfend:end)]);
                    else
                        fieldgain=nan;
                    end
                else
                    fieldgain=nan;
                end
            catch
                fieldgain=nan;
            end
            PFgains(cumct)=fieldgain; cumct=cumct+1;
        end
    end
    
    subplot(1,4,3);
    boxScatterplot((PFgains(celltypes>0)),celltypes(celltypes>0),'xLabels',...
        {'Odor Unresponsive','Odor Responsive','Odor Selective'},'position',[]);
    ylabel('PFgain (rate)');  xtickangle(340); box off;    
    [p,tbl,stats]=kruskalwallis(PFgains(celltypes>0),celltypes(celltypes>0),'off');

    fprintf('anovan f(%d)=%.2f, p=%.4f \n',tbl{2,3},tbl{2,5},tbl{2,6});
    [p,~,stats]=ranksum(PFgains(celltypes==1),PFgains(celltypes>1));
    fprintf('%s PFgains: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p);
    title(sprintf('%s PFgains: Ranksum z=%.2f, p=%.4f \n',region{i},stats.zval,p));
    
    subplot(1,4,4);
    % what about overall rate for odor selective/odor responsive? (0 not 1 is)
    allPFs=cellfun(@(a) any(a(1:2)), {cellPool.FiresDuringRun});
    celltypes2=double(allPFs & cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}));
    
    %celltypes2(allPFs & cellfun(@(a) a(3)==1, {cellPool.OdorSelective}))=2;
    
    cellrates=cellfun(@(a) a(1), {cellPool.meanrate}); % overall mean rate
    boxScatterplot((cellrates),celltypes2+1,'xLabels',...
        {'Odor Unresponsive','Odor Responsive','Odor Selective'},'position',[]);
    ylabel('Mean firing rate');  xtickangle(340); box off;
    
    %[p,tbl,stats]=kruskalwallis(cellrates,celltypes2','off');
    %fprintf('PF rate: anovan f(%d)=%.2f, p=%.4f \n',tbl{2,3},tbl{2,5},tbl{2,6});
    [p,~,stats]=ranksum(cellrates(celltypes2==0),cellrates(celltypes2>0));
    fprintf('%s MeanRate: Ranksum z=%.2f, p=%.2e \n',region{i},stats.zval,p);
    title(sprintf('%s MeanRate: Ranksum z=%.2f, p=%.2e \n',region{i},stats.zval,p));
        
    fprintf('\n');
    % sgtitle doesnt really make sense here, gotta add numbers in xlabels
    %sgtitle(sprintf('Place Fields of %s cells from %s',type, region{i}));
end
%% where do place field peaks fall for odor selective units?
figure;
mycolors=lines(3);
for i=1:length(region)
    % pull pyrs from this brain region that are ODOR SELECTIVE
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) a(3)==1, {SuperUnits.OdorSelective});
    
    cellPool=SuperUnits(cellfilt);
    % now I think i can probably take each unit and replicate it twice for out
    % l and out r and then cull those without place fields
    % first gather the place fields into a super long vector
    allPFs=cell2mat(cellfun(@(a) a.PFmaxpos(1:2), {cellPool.FieldProps},'UniformOutput',false)');
    
    % i think i want to see where the place fields fall for the same run
    % and the opposing run
    %[b,a]=histcounts(allPFs(:,1),0:5:100);
    %plot(mean([a(2:end);a(1:end-1)]),SmoothMat2(b,[10 10],1)./sum(~isnan(allPFs(:,1))),'Color',mycolors( i,:));
    %hold on;
    %[b,a]=histcounts(allPFs(:,2),0:5:100);
    %plot(mean([a(2:end);a(1:end-1)]),SmoothMat2(b,[10 10],1)./sum(~isnan(allPFs(:,1))),'Color',mycolors(i,:)*.6);
    
    [b,a]=histcounts(linearize(allPFs),0:5:100);
    plot(mean([a(2:end);a(1:end-1)]),SmoothMat2(b,[10 10],1)./sum(~isnan(allPFs(:,1))),'Color',mycolors( i,:));
    hold on;
end
mykids=get(gca,'Children');
legend(region);
title(sprintf('Distribution of Field Peaks \n of Odor Responsive Units'));
xlabel('Track Position (% of total distance)');
ylabel(sprintf('Relative Prevalence of \n Field Centers'));

%% what about mean rate along track for odor responsive cells?
figure;
mycolors=lines(3);
for i=1:length(region)
    fullcurve=[];
    % pull pyrs from this brain region that are ODOR RESPONSIVE
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) a(3)==1, {SuperUnits.OdorSelective});
    
    cellPool=SuperUnits(cellfilt);
    % now I think i can probably take each unit and replicate it twice for out
    % l and out r and then cull those without place fields
    % first gather the place fields into a super long vector
    allPFs=cellfun(@(a) nanmean(cell2mat(a(1:2)')), {cellPool.LinPlaceFields},'UniformOutput',false);
    allrates=cell2mat(allPFs(cellfun(@(a) length(a)==99, allPFs))');
    
    % i think i want to see where the place fields fall for the same run
    % and the opposing run
    meanrate=nanmean(allrates(:,1:98)); semrate=SEM(allrates(:,1:98),1);
    plot(meanrate,'Color',mycolors(i,:)); hold on;
    patch([1:length(meanrate) length(meanrate):-1:1],[meanrate+semrate ...
        fliplr(meanrate-semrate)],mycolors(i,:),'FaceAlpha',.4,'EdgeColor','none');
end

mykids=get(gca,'Children');
legend(mykids([4 2]),region);
title(sprintf('Mean Ensemble rate \n of Odor Responsive Units'));
xlabel('Track Position (% of total distance)');
ylabel(sprintf('Ensemble Firing Rate (hz)'));

%%
% do cells that respond to odors respond to place differently?
%{
1. do odor responsive cells have more or less path equivalence across the
two outbound runs, and do they prefer one over the other?  basically xcorr
between runs, and difference in peak rate or mean rate for outbound runs

2. do odor selective cells have a similar/opposite preference during
outbound runs? -mean rate on runs, peak rate on runs, information score
along runs?

%}

colors=lines(2);

% first odor responsive cells
for i=1:length(region)
    % filter by region and celltype
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) (length(a)>1), {SuperUnits.LinPlaceFields});
    cellPool=SuperUnits(cellfilt);
    
    % second filter, has to have rate curves on both outbounds
    cellPool=cellPool(cellfun(@(a)  ~isempty(a{1}) && ~isempty(a{2}), {cellPool.LinPlaceFields}));
    % kill those with too few trajectories

    % lets categorize the cells so we can plot their continuous dims
    % selectivity is -1 0 or 1, y/n and preferred dir
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) diff(a(1:2))>0, {cellPool.OdorMeans})-.5)*2);
    % objres is 1 rate 2 change 3 is dprime 4 is signrankp
    % -1 if downreg, 0 if no reg, 1 if upreg
    objRes=cellfun(@(a) min(a(:,4))<.05,{cellPool.OdorResponsive}).*...
        ((cellfun(@(a) nanmean(a(:,2))>0,{cellPool.OdorResponsive})-.5)*2);
    
    % 0 for no object resp, 1 for resp, and 2 for sel, you can gather direction
    celltypes=double(objRes~=0 & objSel==0);
    celltypes(objSel~=0)=2;
    
    % 1. does object responsiveness relate to corelation between trajectories?
    % first gather the object responsiveness (signed log10)
    % signed effect size
    objscores=cellfun(@(a)  mean(a(:,3))*((mean(a(:,2))>0)-.5)*2, {cellPool.OdorResponsive});
    
    % gather the trajectory correlations too ** some have empty cells
    trajcorrs=cellfun(@(a) corr(a{1}',a{2}','Rows','Pairwise'),{cellPool.LinPlaceFields});
    
    figure; 
    scatter(objscores(celltypes==0),trajcorrs(celltypes==0),12,'filled');
    hold on;
    scatter(objscores(celltypes>0),trajcorrs(celltypes>0),12,'filled');
    xlabel('object responsiveness score (chg from baseline)'); ylabel('outbound pf correlation');
    legend('obj nonresponsive','obj responsive');
    title(sprintf('%s %s cells',region{i}, type)); box off;
    % first filter based on has a trajectory, then run corr on trajectories
    
    
    % Does odor selectivity match trajectory selectivity (mean rates across
    % runs)
    % this is signed dprime effect size
    odordiffs=cellfun(@(a) a(4), {cellPool.OdorSelective}).*...
        cellfun(@(a) (((a(1)-a(2))>0)-.5)*2, {cellPool.OdorMeans});
    % this is dprime effect size betweeen mean rates across runs
    trajdiffs=cellfun(@(a) dprime(a{1},a{2}).*((mean(a{1})>mean(a{2}))-.5)*2, {cellPool.RunRates});
    
    figure; 
    scatter(odordiffs(celltypes<2),trajdiffs(celltypes<2),12,'filled');
    hold on;
    scatter(odordiffs(celltypes==2),trajdiffs(celltypes==2),12,'filled');
    xlabel('d'' obj selectivity'); ylabel('d'' outbout pf mean rate');
    legend('obj nonselective','obj selective');
    title(sprintf('%s %s cells',region{i}, type));
    
    
    
    
    % Does odor selectivity match diff in peak rates fo each run
    % signed d' of odor coding
    odordiffs=cellfun(@(a) a(4), {cellPool.OdorSelective}).*...
        cellfun(@(a) (((a(1)-a(2))>0)-.5)*2, {cellPool.OdorMeans});

    
    % this is just a raw difference in max
    trajdiffs=cellfun(@(a) (max(a{1})-max(a{2}))/...
        (max(a{1})+max(a{2})), {cellPool.LinPlaceFields});
    
    
    
    figure; 
    %sc=scatter(odordiffs(celltypes<2),trajdiffs(celltypes<2),12,'filled'); hold on; 
    %sc(2)=scatter(odordiffs(celltypes==2),trajdiffs(celltypes==2),12,'filled');
    
    sc=scatter(odordiffs,trajdiffs,12,'filled');
    xlabel('Odor Selectivity (d'')'); ylabel('Outbound Trajectory Selectivity (grand mean diff/sum)');
    
    hold on;
    plot([-1 1],[0 0],'k'); plot([0 0],[-1 1],'k');
    
    xlim([-1 1]); ylim([-1 1]);
    [rho,p]=corr(odordiffs',trajdiffs','rows','complete');
    [rho2,p2]=corr(odordiffs(celltypes==2)',trajdiffs(celltypes==2)','rows','complete');
    fprintf('for full %s pop r^2=%.2f, p=%.4f, for selective cells, r^2=%.2f, p=%.4f \n',...
        region{i},rho,p,rho2,p2);
    
    [r,m,b] = regression(odordiffs',trajdiffs','one');
    xlims=get(gca,'XLim'); plot(xlims, xlims.*m+b,'color',colors(i,:));
    title(sprintf('%s %s cells r^2=%.2f, p=%.4f',region{i}, type, rho, p));
    legend(sc,'All units','Regression Line'); box off;
    % howabout whether the cell had a place field ont he same or opposing
    % trajectory
    
    
    %somethiing like none, same, opposite both
    % we need object selectivity, and the direction
    % thats objSel,  now we need pf presence, so it will have to be something
    % like -1, 0, 1, or 2, and then if it has one, we multiply it by the diff
    pfnums=cellfun(@(a) sum(a(1:2)), {cellPool.PFexist});
    pfdir=cellfun(@(a) diff(a(1:2)), {cellPool.PFexist});
    pfnums(pfdir~=0)=pfnums(pfdir~=0).*pfdir(pfdir~=0); % sign all those who are 1
    % 53 left, 300 no pfs, 73 both, and 57 right
    % now aggregate,
    
    bary=[nanmean(pfnums(objSel~=0)==0),...                  % sum no pf
        nanmean(pfnums(objSel~=0)==objSel(objSel~=0)),... % sum same side
        nanmean(pfnums(objSel~=0)==-objSel(objSel~=0)),... % sum opposite side
        nanmean(pfnums(objSel~=0)==2)];                   % sum both sides
    figure;
    b=bar(bary);
    set(gca,'XTickLabel',{sprintf('No Place fields'),...
        sprintf('PF in matching route'),...
        sprintf('Opposite route'),...
        sprintf('Both routes')});
    %legend('Object Responsive','Object Selective');
    ylabel('% of units');
    xtips = b(1).XEndPoints; ytips = b(1).YEndPoints;
    labels = {sum(pfnums(objSel~=0)==0),...
        sum(pfnums(objSel~=0)==objSel(objSel~=0)),...
        sum(pfnums(objSel~=0)==-objSel(objSel~=0)),...
        sum(pfnums(objSel~=0)==2)};
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    title(sprintf('%s pf odds %d same vs opposite',region{i},1));
    xtickangle(340); box off;
    
    % probably a binomial against 50%?
    pfcounts=cell2mat(labels);
    binoprob=binocdf(max(pfcounts(:,2:3)),sum(pfcounts(2:3)),0.5);
    fprintf('Bino test on same side vs opposite side pf %.4f \n',1-binoprob);
    
    
    % howabout a boxscatter with x as left odor preferring and right odor
    % preferring, and y trajectory preference
    
    % add 1 to these to get our vals
    objSel=cellfun(@(a) a(3)==1,{cellPool.OdorSelective}).*...
        ((cellfun(@(a) diff(a(1:2))>0, {cellPool.OdorMeans})-.5)*2);
    
    % turn right preferring into 2
    objSel(objSel==-1)=2;
    % diff over mean of average firing rates for each run
    trajdiffs=cellfun(@(a) (nanmean(a{1})-nanmean(a{2}))/...
        (nanmean(a{1})+nanmean(a{2})), {cellPool.RunRates});
    nanout=isnan(objSel) |isnan(trajdiffs);
    boxScatterplot(trajdiffs(objSel>0 & ~nanout)',objSel(objSel>0 & ~nanout),'xLabels',{'Left Preferring','Right Preferring'});
    ylabel({'Trajectory Selectivity'; 'Left Preferring   Right Preferring'});
    xlabel('Odor Period Preference');
    [p,h,stats]=ranksum(trajdiffs(objSel==1 & ~nanout),trajdiffs(objSel==2 & ~nanout));
    title(sprintf('%s ranksum test \n L mean=%.2f R mean %.2f \n z=%.2f, p=%.2e',...
        region{i},nanmean(trajdiffs(objSel==1)),nanmean(trajdiffs(objSel==2)),stats.zval,p));
    set(gca,'YLim',[-1 1]);
    
end
%% the final odor vs route figure by field tercile

% first odor responsive cells
for i=1:length(region)
    % filter by region and celltype
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) (length(a{1})>1) & length(a{2})>1, {SuperUnits.RunRates});
    cellPool=SuperUnits(cellfilt);
    
    % second filter, has to have rate curves on both outbounds
    cellPool=cellPool(cellfun(@(a)  max(a{1})>0 && max(a{2})>0, {cellPool.LinPlaceFields}));
    
    
    % now find the position of the peak
    peakpos=cellfun(@(a) a.PFmaxpos(find(a.PFmax(1:2)==max(a.PFmax(1:2)),1,'first')),...
        {cellPool.FieldProps}, 'UniformOutput', false);

    
    % 1. does object responsiveness relate to corelation between trajectories?
    % first gather the object responsiveness (signed log10)
    
    % ODOR (dprime or diff over sum of means)
    % diff over sum
    % signed d' of odor coding (dprime, sign by which mean is larger)
    %odordiffs=cellfun(@(a) a(4), {cellPool.OdorSelective}).*cellfun(@(a) (((a(1)-a(2))>0)-.5)*2, {cellPool.OdorMeans});
    myxlabel='d'' Odor selectivity';
    odordiffs=cellfun(@(a) (a(1)-a(2))/sum(a), {cellPool.OdorMeans});
    myxlabel='Odor Selectivity Index';
    
    
    % TRAJECTORY (dprime sum over run, diff/sum peak, or diff/sum mean across
    % space or runs
    % gather the trajectory correlations too ** some have empty cells
    % correlation?
    % trajcorrs=cellfun(@(a) corr(a{1}',a{2}','Rows','Pairwise'),{cellPool.LinPlaceFields});
    % this is just diff over sum max across space
    %trajdiffs=cellfun(@(a) (max(a{1})-max(a{2}))/...
    %    (max(a{1})+max(a{2})), {cellPool.LinPlaceFields});
    
    % this is the difference in RunRates
    %trajdiffs=cellfun(@(a) a(1), {cellPool.SplitterScore});
    
    % diff over sum over mean rate across runs
    trajdiffs=cellfun(@(a) (nanmean(a{1})-nanmean(a{2}))/...
        (nanmean(a{1})+nanmean(a{2})), {cellPool.RunRates});
    
    %trajdiffs=cellfun(@(a) dprime(a{1},a{2}), {cellPool.RunRates});
    % diff over sum over peak rates on trajs
    %trajdiffs=cellfun(@(a) (max(a{1})-max(a{2}))/...
    %    (max(a{1})+max(a{2})), {cellPool.LinPlaceFields});
    myylabel='Outbound Trajectory Selectivity';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % have to grab the position of the peak, and use the pos of the highest
    % peak
    peakpos(cellfun(@(a) isempty(a), peakpos))={nan};
    peakpos=cell2mat(peakpos);
    figure;
    %posinds=[1 20; 21 40; 41 60; 61 80; 81 100]; mylegend={'tercile 1',
    % ...'tercile 2','tercile 3','tercile 4','tercile 5'};
    %posinds=[1 25; 26 50; 51 75; 76 100];
    %posinds=[1 33; 34 66; 67 100];
    posinds=[1 20; 20 101]; mylegend={'First tercile','Last four'};
    %posinds=[1 100]; mylegend={'First tercile','Last four'};0.1736
    
    
    %colors=hsv(size(posinds,1));
    colors=lines(size(posinds,1));
    %colors=winter(size(posinds,1));
    %colors=[.7 .2 .1; 0 0 0];
    
    sc=[];
    for posi=1:size(posinds,1)
        subplot(1,size(posinds,1)+1,posi);
        sc(posi)=scatter(odordiffs(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2)),...
            trajdiffs(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2)),...
            12,colors(posi,:),'filled');
        hold on;
        [r,pval]=corr(odordiffs(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2))',...
            trajdiffs(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2))','rows','complete');
        fprintf('%s bin %d, r=%.2f p=%.4f \n', region{i}, posi,r,pval);
        
        
        xlabel(myxlabel); ylabel(myylabel);
        title(sprintf('%s %s Peaks at %s,  %d cells r=%.2f p=%.4f',...
            region{i}, type,num2str(posinds(posi,:)),...
            sum(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2) & ~isnan(odordiffs) & ~isnan(trajdiffs)),...
            r,pval));
        
        hold on; 
        myx=abs(max(get(gca,'XLim'))); plot([-myx myx],[0 0],'k'); 
        myy=abs(max(get(gca,'YLim'))); plot([0 0],[-myy myy],'k');
        
        [r,m,b] = regression(odordiffs(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2))',...
            trajdiffs(peakpos>=posinds(posi,1) & peakpos<posinds(posi,2))','one');
        plot([-myx myx], [-myx myx].*m+b,'color',colors(posi,:));
        %legend(sc,mylegend); box off;
        xlim([-myx myx]);
        ylim([-myy myy]);
    end
    
    subplot(1,size(posinds,1)+1,size(posinds,1)+1);
    scatter(odordiffs,trajdiffs,12,[.6 .6 .6],'filled'); hold on;
    
    myx=abs(max(get(gca,'XLim'))); plot([-myx myx],[0 0],'k');
    myy=abs(max(get(gca,'YLim'))); plot([0 0],[-myy myy],'k');
    
    [r,m,b] = regression(odordiffs(peakpos>=posinds(1) & peakpos<posinds(end,end))',...
        trajdiffs(peakpos>=posinds(1) & peakpos<posinds(end,end))','one');
    plot([-myx myx], [-myx myx].*m+b,'color',[.5 .5 .5]);
    %legend(sc,mylegend); box off;
    xlim([-myx myx]);
    ylim([-myy myy]);
    [r,pval]=corr(odordiffs',trajdiffs','rows','complete');
    fprintf('%s all bins, %d cells, r=%.2f p=%.4f \n', region{i},r,pval);
     
    title(sprintf('%s %s tot:%d cells r=%.2f p=%.4f',region{i},type,...
         sum(peakpos>=posinds(1) & peakpos<posinds(end,end) & ~isnan(odordiffs) & ~isnan(trajdiffs)),...
         r,pval));
     xlabel(myxlabel); ylabel(myylabel);
     set(gcf,'Position',[2 558 1533 420]);

end

fprintf('\n');
%%
% shantanu had a question, basically the question is whether route
% selective cells, or cells whose place fields significantly differ across
% the routes behave differently than cells who arent...
% the baseline here i think is cells that dont?

% first gather the cells whose pf maps differ
% those are the cells whose SplitterScore is nonzero and p is less than 0


% the other problem is that it might be better to filter splitter cells to
% cells who fire on the center stem and arent head direction cells, or
% distal arm cells.. whatever, lets look for cells that have pfs but arent
% splitters and then compare those to splitters
region={'CA1','PFC'};
type='pyr';
for i=1:length(region)
    
    figure;
    % filter by region and celltype and has to fire during at least one run AND
    % has to fire during odor delivery
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) any(a),{SuperUnits.FiresDuringRun});
    
    % of the splitter cells, whats their object selectivity?
    cellPool=SuperUnits(cellfilt);
    Splitters=cellfun(@(a) a(1)>0 & a(2)<.01, {cellPool.SplitterScore}); % pval splitter
    
    % now compare rates of odor responsiveness and odor selectivity of the
    % splitters to the non splitters
    OdorSelectivity=cellfun(@(a) a(1), {cellPool.OdorSelective}); % diff/sum
    subplot(1,5,1);
    boxScatterplot(OdorSelectivity,Splitters+1,'Position',[],'xLabels',{'nonSplitter','Splitter'});
    ylabel('Odor Selectivity Score'); ylim([0 1])
    [p,~,stats]=ranksum(OdorSelectivity(Splitters==0),OdorSelectivity(Splitters==1));
    title(sprintf('%s odor selectivity \n Z(%d)=%.2f, p=%.3f)',region{i},stats.ranksum,stats.zval,p));
    
    
    % how about odor rate change
    subplot(1,5,2);
    OdorResponses=cell2mat(cellfun(@(a) a(:,2), {cellPool.OdorResponsive}, 'UniformOutput', false)'); % delta rate
    boxScatterplot(OdorResponses,linearize(repmat(Splitters'+1,2,1)),'Position',[],'xLabels',{'nonSplitter','Splitter'});
    ylabel('Odor Induced Rate Change');
    
    [p,~,stats]=ranksum(OdorResponses(linearize(repmat(Splitters',2,1))==0),...
        OdorResponses(linearize(repmat(Splitters',2,1))==1));
    title(sprintf('%s odor responsiveness \n Z(%d)=%.2f, p=%.3f)',region{i},stats.ranksum,stats.zval,p));
    
    
    % likelihood of being significantly odor responsive, and significantly odor
    % selective
    isOdorResponsive=cellfun(@(a) any(a(:,4)<.05), {cellPool.OdorResponsive});
    isOdorSelective=cellfun(@(a) a(3)==1, {cellPool.OdorSelective});
    
    stackedy=[nanmean(isOdorResponsive(Splitters==0) & ~isOdorSelective(Splitters==0)),...
        nanmean(isOdorSelective(Splitters==0));...
        nanmean(isOdorResponsive(Splitters==1) & ~isOdorSelective(Splitters==1)),...
        nanmean(isOdorSelective(Splitters==1))]*100;
    
    
    subplot(1,5,3);
    b=bar(stackedy,'stacked');
    set(gca,'XTickLabel',{'nonSplitter','Splitter'});
    legend('Odor Responsive','Odor Selective');
    ylabel('% of units');
    
    xtips = b(2).XEndPoints; ytips = b(1).YEndPoints;  ytips2 = b(2).YEndPoints;
    labels = {sum(isOdorResponsive(Splitters==0) & ~isOdorSelective(Splitters==0)),...
        sum(isOdorResponsive(Splitters==1) & ~isOdorSelective(Splitters==1))};
    
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    
    labels2 =  {sum(isOdorSelective(Splitters==0));...
        sum(isOdorSelective(Splitters==1))};
    
    text(xtips,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    
    chistats=[sum(~isOdorResponsive(Splitters==0)),...
        sum(isOdorResponsive(Splitters==0) & ~isOdorSelective(Splitters==0)),...
        sum(isOdorSelective(Splitters==0));...
        sum(~isOdorResponsive(Splitters==1)),...
        sum(isOdorResponsive(Splitters==1) & ~isOdorSelective(Splitters==1)),...
        sum(isOdorSelective(Splitters==1))];
    [pval,chi2stat,df]=chi2indep(chistats');
    title(sprintf('Selectivity rates, \n Chi(%d)=%.2f, p=%.3f',df,chi2stat,pval));
    
    [pval,chi2stat,df]=chi2indep(chistats(:,[1 2])');
    fprintf('%s Odor Responsive rates Splitters, Chi(%d)=%.2f, p=%.3f \n',region{i},df,chi2stat,pval);
    
    [pval,chi2stat,df]=chi2indep(chistats(:,[1 3])');
    fprintf('%s Odor Selective rates Splitters, Chi(%d)=%.2f, p=%.3f \n',region{i},df,chi2stat,pval);
    
    subplot(1,5,4);
    stackedy=[nanmean(isOdorResponsive(Splitters==0));...
        nanmean(isOdorResponsive(Splitters==1))]*100;
    b=bar(stackedy,'stacked');
    set(gca,'XTickLabel',{'nonSplitter','Splitter'});
    %legend('Odor Responsive','Odor Selective');
    ylabel('% Odor Responsive Units');
    
    xtips = b.XEndPoints; ytips = b.YEndPoints; 
    labels = {sum(isOdorResponsive(Splitters==0)),...
        sum(isOdorResponsive(Splitters==1))};
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    
    binostat=[sum(~isOdorResponsive(Splitters==0)),...
        sum(isOdorResponsive(Splitters==0));...
        sum(~isOdorResponsive(Splitters==1)),...
        sum(isOdorResponsive(Splitters==1))];
    pval=1-binocdf(sum(isOdorResponsive(Splitters==1)),sum(Splitters==1),...
        nanmean(isOdorResponsive(Splitters==0)));
    title(sprintf('Odor responsiveness, \n bino(%d / %d) null=%.f%%, p=%.4f',...
        sum(isOdorResponsive(Splitters==1)),sum(Splitters==1),...
        nanmean(isOdorResponsive(Splitters==0))*100,pval));
    
    
    subplot(1,5,5);
    myrates=cellfun(@(a) a, {cellPool.meanrate});
    boxScatterplot(myrates,Splitters'+1,'Position',[]);
    title(sprintf('%s Mean rate', region{i}));
    set(gca,'XTickLabel',{'nonSplitter','Splitter'});
    [p,~,stats]=ranksum(myrates(Splitters),myrates(~Splitters));
    title(sprintf('%s odor rates \n Z(%d)=%.2f, p=%.3f)',region{i},stats.ranksum,stats.zval,p));
    ylabel('Mean Firing Rate (hz)');
    set(get(gcf,'Children'),'Box','off')

    
end
% so this is kindof a nonstarter, mostly because these cells all split
% based on this metric... maybe i'll need to rework the stats on this, i
% kindof hastily did the splitter score code


