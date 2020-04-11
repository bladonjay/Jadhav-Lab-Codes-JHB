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
 
 
chiraw=[sum(sum(PFcounts(objSel==0 & objRes==0,:),2)==0)...
    sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)==0)...
    sum(sum(PFcounts(objSel~=0,:),2)==0);...
    sum(sum(PFcounts(objSel==0 & objRes==0,:),2)==1)...
    sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)==1)...
    sum(sum(PFcounts(objSel~=0,:),2)==1);...
    sum(sum(PFcounts(objSel==0 & objRes==0,:),2)>1)...
    sum(sum(PFcounts(objSel==0 & objRes~=0,:),2)>1)...
    sum(sum(PFcounts(objSel~=0,:),2)>1)];
[a,b,c]=chi2indep(chiraw);
fprintf('%s Chi2 test, # routs with spikes across odor repsonses chi2(%d)= %.2f p=%.2e \n',...
    region{r},c,b,a)
end
%% okay now to plot odor period rates vs run rates

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
fprintf('%s pyr odor rate to run rate corr=%.2f p=%.2e, HC in corr=%.2f p=%.2e\n',region{i},a,b,a2,b2)

sp(i)=subplot(2,2,i);
scatter(odormeans(regionfilt & PYRfilt & nanfilt),runmeans(regionfilt & PYRfilt & nanfilt),4,mycolors(1,:),'filled');
xlabel('Rate during odorpoke'); ylabel('Rate during run');  set(gca,'Xscale','log','Yscale','log');
title(sprintf('%s r2=%.2f p=%.2e',region{i},a,b)); legend('pyr');
si(i)=subplot(2,2,i+2);
scatter(odormeans(regionfilt & INfilt & nanfilt),runmeans(regionfilt & INfilt & nanfilt),4,mycolors(2,:),'filled'); 
xlabel('Rate during odorpoke'); ylabel('Rate during run');set(gca,'Xscale','log','Yscale','log');
title(sprintf('%s r2=%.2f p=%.2e',region{i},a2,b2)); legend('IN');


end  
linkaxes(sp);
linkaxes(si);


 %%
 % now filter out the cells who dont spike during runs
 
 % I think these cells will be called 'task relevant firing' e.g. fired on
 % tyhe run, or fired during the odor sampling epoch (at least N spikes)
 for i=1:length(region)
 
cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
    cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
    cellfun(@(a) any(a),{SuperUnits.FiresDuringRun});

% just take the center tercile?
usetercile=0;
if usetercile
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type});
    
    overallrate=log2(cellfun(@(a) a, {SuperUnits(cellfilt).meanrate}));
    figure; hist(overallrate);
    lowfilt=prctile(overallrate,33); highfilt=prctile(overallrate,66);
    % just take the center tercile?
    cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
        cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
        cellfun(@(a) log2(a)>lowfilt & log2(a)<highfilt, {SuperUnits.meanrate});
end

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
cellclasses=double(objSel==0 & objRes==0); cellclasses(objSel==0 & objRes~=0)=2; cellclasses(objSel~=0)=3;

% ok %% of cells with 1 and %% of cells with more than 1 pf that arent obj
PFcounts=cell2mat(cellfun(@(a) a(1:2), {cellPool.PFexist},'UniformOutput',false)');

bary=[nanmean(sum(PFcounts(cellclasses==1,:),2)==1)...          % nonselective
    nanmean(sum(PFcounts(cellclasses==2,:),2)==1)...            % odor responsive
    nanmean(sum(PFcounts(cellclasses==3,:),2)==1);...           % odor selective
    nanmean(sum(PFcounts(cellclasses==1,:),2)>1)...
    nanmean(sum(PFcounts(cellclasses==2,:),2)>1)...
    nanmean(sum(PFcounts(cellclasses==3,:),2)>1)]*100;
figure;
b=bar(bary','stacked');
set(gca,'XTickLabel',{sprintf('Non Odor Responsive (%d)',sum(cellclasses==1)),...
     sprintf('Odor Responsive (%d)', sum(cellclasses==2)),...
     sprintf('Odor Selective (%d)', sum(cellclasses==3))});
legend('PF on one outbound rounte','PF on both outbounds');
ylabel('% of units with fields');
xtips = b(2).XEndPoints; ytips = b(1).YEndPoints;  ytips2 = b(2).YEndPoints; 
labels = {sum(sum(PFcounts(cellclasses==1,:),2)==1) sum(sum(PFcounts(cellclasses==2,:),2)==1) sum(sum(PFcounts(cellclasses==3,:),2)==1)};
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
labels2 = {sum(sum(PFcounts(cellclasses==1,:),2)>1) sum(sum(PFcounts(cellclasses==2,:),2)>1) sum(sum(PFcounts(cellclasses==3,:),2)>1)};
text(xtips,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
ylim([0 100]);
sgtitle(sprintf('Place Fields of %s cells from %s, out of a pool of %d active', type, region{i}, length(cellPool)));
xtickangle(340); box off;


% raw numbers and their perceotages
fprintf('%s numbers: Nonresponsive tot %d: one %d (%.f): two %d (%.f) \n',...
   region{i},sum(objSel==0 & objRes==0), labels{1}, bary(1,1), labels2{1},bary(2,1));
 fprintf('%s numbers: Responsive tot %d: one %d (%.f): two %d (%.f) \n',...
   region{i},sum(objSel==0 & objRes~=0), labels{2}, bary(1,2), labels2{2},bary(2,2));
 fprintf('%s numbers: Selective tot %d: one %d (%.f): two %d (%.f) \n',...
   region{i},sum(objSel~=0), labels{3}, bary(1,3), labels2{3}, bary(2,3));
 


% chi square test for independence
chiraw=[sum(sum(PFcounts(cellclasses==1,:),2)==0)...        % nonselective
    sum(sum(PFcounts(cellclasses==2,:),2)==0)...            % odor responsive
    sum(sum(PFcounts(cellclasses==3,:),2)==0);...           % odor selective
    sum(sum(PFcounts(cellclasses==1,:),2)==1)...            % nonselective
    sum(sum(PFcounts(cellclasses==2,:),2)==1)...            % odor responsive
    sum(sum(PFcounts(cellclasses==3,:),2)==1);...           % odor selective
    sum(sum(PFcounts(cellclasses==1,:),2)>1)...
    sum(sum(PFcounts(cellclasses==2,:),2)>1)...
    sum(sum(PFcounts(cellclasses==3,:),2)>1)];
[a,b,c]=chi2indep(chiraw);

% the chi2 tests...
fprintf('%s Chi2 on pf rates across object responsiveness is chi2(%d)=%.2f, p%.2e \n',...
    region{i},c,b,a);
[a,b,c]=chi2indep(chiraw(:,1:2));
fprintf('%s posthoc unresp v resp chi2(%d)=%.2f, bonf. p=%.2e \n',...
    region{i},c,b,a*3);
[a,b,c]=chi2indep(chiraw(:,[1 3]));
fprintf('%s posthoc unresp v sel chi2(%d)=%.2f, bonf. p=%.2e \n',...
    region{i},c,b,a*3);
[a,b,c]=chi2indep(chiraw(:,2:3));
fprintf('%s posthoc resp v sel chi2(%d)=%.2f, bonf. p=%.2e \n',...
    region{i},c,b,a*3);

% now run binomial test on pf prevalence
[y]=binocdf(chiraw(2,2)+chiraw(3,2),sum(chiraw(:,2)),(chiraw(2,1)+chiraw(3,1))/sum(chiraw(:,1)));
fprintf('binomial test on p(odorresponsive) cells  pval=%.2e\n',1-y);
[y]=binocdf(chiraw(2,3)+chiraw(3,3),sum(chiraw(:,3)),(chiraw(2,1)+chiraw(3,1))/sum(chiraw(:,1)));
fprintf('binomial test on p(odorselective) cells  pval=%.2e\n',1-y)

 end
 %%
 for i=1:length(region)
ratecontrol=1;
if ratecontrol
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
    cellclasses=double(objSel==0 & objRes==0); cellclasses(objSel==0 & objRes~=0)=2; cellclasses(objSel~=0)=3;
    PFcounts=cell2mat(cellfun(@(a) a(1:2), {cellPool.PFexist},'UniformOutput',false)');

   
    % this could be the overall firing rate...
    overallrate=(cellfun(@(a)  a, {cellPool.meanrate}));
    
    meansy=[nanmean(overallrate(cellclasses==1)) nanmean(overallrate(cellclasses==2))...
        nanmean(overallrate(cellclasses==3))];
    SEMsy=[SEM(overallrate(cellclasses==1)) SEM(overallrate(cellclasses==2))...
        SEM(overallrate(cellclasses==3))];
    figure; bar(1:3,meansy); hold on; errorbar(1:3, meansy, SEMsy,'k.');
    set(gca,'XTickLabel',{'Task Responsive','Odor Responsive','Odor Selective'});
    xtickangle(-10); ylabel('Log_2 Mean Rate');
    sgtitle(region{i});
    
    meansy=[nanmean(overallrate(sum(PFcounts,2)==0)) nanmean(overallrate(sum(PFcounts,2)==1))...
        nanmean(overallrate(sum(PFcounts,2)==2))];
    SEMsy=[SEM(overallrate(sum(PFcounts,2)==0)) SEM(overallrate(sum(PFcounts,2)==1))...
        SEM(overallrate(sum(PFcounts,2)==2))];
    figure; bar(1:3,meansy); hold on; errorbar(1:3, meansy, SEMsy,'k.');
    set(gca,'XTickLabel',{'No PFs','One PF','2+ PFs'});
    xtickangle(340); ylabel(' Mean Rate');
    sgtitle(region{i});
    % stats are friedman test if you want them...
end
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

chiraw=[sum(objRes(PFcounts==0)==0)...
    sum(objRes(PFcounts==1)==0)...
    sum(objRes(PFcounts>1)==0);...
    sum(objRes(PFcounts==0)~=0 & objSel(PFcounts==0)==0)...  % %% of 0 pf units who are responsive        
    sum(objRes(PFcounts==1)~=0 & objSel(PFcounts==1)==0)...   % %% of 1 pf units who are responsive
    sum(objRes(PFcounts>1)~=0 & objSel(PFcounts>1)==0);...  % %% of 2 pf units who are responsive               
    sum(objSel(PFcounts==0)~=0)...  
    sum(objSel(PFcounts==1)~=0)...    % %% of 1 pf units who are selective
    sum(objSel(PFcounts>1)~=0)];

[a,b,c]=chi2indep(chiraw);
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



% this plots ECDFS

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
    objSel=objSel';
    % -1 if downreg, 0 if no reg, 1 if upreg
    objRes=repmat(cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}),2,1); %.*...
    %((cellfun(@(a) a(2)>0,{PFtest.OdorResponsive})-.5)*2),2,1);
    objRes=objRes';
    
    % how many cells w pfs are odor selective or responsive...
    celltypes=double(allPFs(:)==1 & objSel(:)==0 & objRes(:)==0);
    celltypes(allPFs(:)==1 & objRes(:)~=0 & objSel(:)==0)=2;
    celltypes(allPFs(:)==1 & objSel(:)~=0)=3;
    
    figure;
    % first is pf peak rate
    PFpeaks=cell2mat(cellfun(@(a) a.PFmax(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    subplot(1,2,1);
    boxScatterplot(PFpeaks(celltypes>0),celltypes(celltypes>0),'xLabels',...
        {'Obj Unresponsive','Obj Responsive','Obj Selective'},'position',[]);
    ylabel('%s Place Field Peak Rate (Hz)');  xtickangle(340); box off;
    [p,tbl,stats]=anovan(PFpeaks(celltypes>0),celltypes(celltypes>0),'display','off');
    [p,tbl,stats]=kruskalwallis(PFpeaks(celltypes>0),celltypes(celltypes>0),'off');
    [c,m] = multcompare(stats,'Display','off');
    fprintf('%s Peak Rate: KruskalWallis Chi2(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    fprintf('%s Nonresponsive vs. responsive diff= %.2f, p=%.4f \n',region{i},c(1,4),c(1,6));
    fprintf('%s Nonresponsive vs. selective diff= %.2f, p=%.4f \n',region{i},c(2,4),c(2,6));
    fprintf('%s responsive vs. responsive diff= %.2f, p=%.4f \n',region{i},c(3,4),c(3,6));
    fprintf('\n');
    
    
    % this does boxscatterplots, mainly cause they make everything look the
    % Place field Width
    subplot(1,2,2);
    PFwidth=cell2mat(cellfun(@(a) a.PFsize(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    boxScatterplot(PFwidth(celltypes>0),celltypes(celltypes>0),'xLabels',...
        {'Obj Unresponsive','Obj Responsive','Obj Selective'},'Position',[]);
    ylabel('Place Field Width (% of maze)');  xtickangle(340); box off;
    sgtitle(sprintf('Place Fields of %s cells from %s',type, region{i}));
    [p,tbl,stats]=anovan(PFwidth(celltypes>0),celltypes(celltypes>0),'display','off');
    [p,tbl,stats]=kruskalwallis(PFwidth(celltypes>0),celltypes(celltypes>0),'off');
    [c,m] = multcompare(stats,'Display','off');
    fprintf('%s PF Width: anovan f(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    fprintf('%s Unresponsive vs responsive diff= %.2f, p=%.4f \n',region{i},c(1,4),c(1,6));
    fprintf('%s Unresponsive vs selective zdiff= %.2f, p=%.4f \n',region{i},c(2,4),c(2,6));
    fprintf('%s responsive vs responsive diff= %.2f, p=%.4f \n',region{i},c(3,4),c(3,6));
    fprintf('\n');
    
    % next is sparsity
    figure;
    subplot(1,4,1);
    PFsparsity=cell2mat(cellfun(@(a) a.sparsity(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    boxScatterplot(log2(PFsparsity(celltypes>0)),celltypes(celltypes>0),'xLabels',...
        {'Obj Unresponsive','Obj Responsive','Obj Selective'},'position',[]);
    ylabel('Log_2 Place Field sparsity');  xtickangle(340); box off;
    [p,tbl,stats]=anovan(PFsparsity(celltypes>0),celltypes(celltypes>0),'display','off');
    [p,tbl,stats]=kruskalwallis(PFsparsity(celltypes>0),celltypes(celltypes>0),'off');
    [c,m] = multcompare(stats,'Display','off');
    fprintf('%s PF Sparsity: KruskalWallis f(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    fprintf('%s Unresponsive vs responsive z= %.2f, p=%.4f \n',region{i},c(1,4),c(1,6));
    fprintf('%s Unresponsive vs selective z= %.2f, p=%.4f \n',region{i},c(2,4),c(2,6));
    fprintf('%s responsive vs responsive z= %.2f, p=%.4f \n',region{i},c(3,4),c(3,6));
    fprintf('\n');
    
    % next is information in bits/spike
    subplot(1,4,2);
    PFinfo=cell2mat(cellfun(@(a) a.info(1:2), {cellPool.FieldProps}, 'UniformOutput', false));
    boxScatterplot(log2(PFinfo(celltypes>0)),celltypes(celltypes>0),'xLabels',...
        {'Obj Unresponsive','Obj Responsive','Obj Selective'},'position',[]);
    ylabel('Log_2 Place Field information (bits/spike)'); xtickangle(340); box off;
    
    [p,tbl,stats]=anovan(PFinfo(celltypes>0),celltypes(celltypes>0),'display','off');
    [p,tbl,stats]=kruskalwallis(PFinfo(celltypes>0),celltypes(celltypes>0),'off');
    [c,m] = multcompare(stats,'Display','off');
    fprintf('%s PF Information: anovan f(%d)=%.2f, p=%.4f \n',region{i},tbl{2,3},tbl{2,5},tbl{2,6});
    fprintf('%s Unresponsive vs responsive z= %.2f, p=%.4f \n',region{i},c(1,4),c(1,6));
    fprintf('%s Unresponsive vs selective z= %.2f, p=%.4f \n',region{i},c(2,4),c(2,6));
    fprintf('%s responsive vs responsive z= %.2f, p=%.4f \n',region{i},c(3,4),c(3,6));
    fprintf('\n');
    
    % next is field gain (this is a hard one...)
    % here though, I am only using cells whose field isnt on the edge...
    cumct=1;
    for j=1:length(cellPool)
        try
            linfield=cellPool(j).LinPlaceFields{1}(1,:);
            [pfmax,pfmaxpos]=max(linfield);
            pfstart=find(linfield(1:pfmaxpos) < pfmax*.25,1,'last');
            if isempty(pfstart), pfstart=1; end
            pfend=pfmaxpos+find(linfield(pfmaxpos:end) < pfmax*.25,1,'first');
            if isempty(pfend), pfend=99; end
            if ~isnan(pfend) && ~isnan(pfstart)
                fieldgain=nanmean(linfield(pfstart+1:pfend-1))-nanmean([linfield(1:pfstart) linfield(pfend:end)]);
            else
                fieldgain=nan;
            end
        catch
            fieldgain=nan;
        end
        PFgains(cumct)=fieldgain; cumct=cumct+1;
        
        try
            linfield=cellPool(j).LinPlaceFields{1}(2,:);
            [pfmax,pfmaxpos]=max(linfield);
            pfstart=find(linfield(1:pfmaxpos) < pfmax*.25,1,'last');
            if isempty(pfstart), pfstart=nan; end
            pfend=pfmaxpos+find(linfield(pfmaxpos:end) < pfmax*.25,1,'first');
            if isempty(pfend), pfend=nan; end
            if ~isnan(pfend) && ~isnan(pfstart)
                fieldgain=nanmean(linfield(pfstart+1:pfend-1))-nanmean([linfield(1:pfstart) linfield(pfend:end)]);
            else
                fieldgain=nan;
            end
        catch
            fieldgain=nan;
        end
        PFgains(cumct)=fieldgain; cumct=cumct+1;
    end
    
    subplot(1,4,3);
    boxScatterplot((PFgains(celltypes>0)),celltypes(celltypes>0),'xLabels',...
        {'Obj Unresponsive','Obj Responsive','Obj Selective'},'position',[]);
    ylabel('PFgain (rate)');  xtickangle(340); box off;
    
    [p,tbl,stats]=anovan(PFgains(celltypes>0),celltypes(celltypes>0),'display','off');
    fprintf('%s PF gain: unresponsive:%.2f +/- %.2f   Resp:%.2f +/- %.2f   Select:%.2f +/- %.2f \n',region{i},...
        nanmean(PFgains(celltypes==1)),SEM(PFgains(celltypes==1)),...
        nanmean(PFgains(celltypes==2)),SEM(PFgains(celltypes==2)),...
        nanmean(PFgains(celltypes==3)),SEM(PFgains(celltypes==3)));
    fprintf('anovan f(%d)=%.2f, p=%.4f \n',tbl{2,3},tbl{2,6},tbl{2,7});
    [p,h,stats]=ranksum(PFgains(celltypes==1),PFgains(celltypes==2),'method','approximate');
    fprintf('Unresponsive vs responsive z= %.2f, p=%.4f \n',stats.zval,p);
    [p,h,stats]=ranksum(PFgains(celltypes==1),PFgains(celltypes==3),'method','approximate');
    fprintf('Unresponsive vs selective z= %.2f, p=%.4f \n',stats.zval,p);
    [p,h,stats]=ranksum(PFgains(celltypes==2),PFgains(celltypes==3),'method','approximate');
    fprintf('responsive vs responsive z= %.2f, p=%.4f \n',stats.zval,p);
    fprintf('\n');
    
    subplot(1,4,4);
    % what about overall rate for odor selective/odor responsive?
    allPFs=cellfun(@(a) any(a(1:2)), {cellPool.FiresDuringRun});
    celltypes2=double(allPFs & cellfun(@(a) any(a(:,4)<.05),{cellPool.OdorResponsive}));
    celltypes2(allPFs & cellfun(@(a) a(3)==1, {cellPool.OdorSelective}))=2;
    
    cellrates=cellfun(@(a) a(1), {cellPool.meanrate}); % overall mean rate
    boxScatterplot(log2(cellrates),celltypes2+1,'xLabels',...
        {'Obj Unresponsive','Obj Responsive','Obj Selective'},'position',[]);
    ylabel('Log_2 cell firing rate');  xtickangle(340); box off;
    
    [p,tbl,stats]=anovan(cellrates,celltypes2','display','off');
    fprintf('%s mean rate unresp mean=%.2f, resp mean=%.2f, & sel mean=%,2f \n',region{i},...
        nanmean(cellrates(celltypes2==0)),nanmean(cellrates(celltypes2==1)),...
        nanmean(cellrates(celltypes2==2)));
    fprintf('PF rate: anovan f(%d)=%.2f, p=%.4f \n',tbl{2,3},tbl{2,6},tbl{2,7});
    [p,h,stats]=ranksum(cellrates(celltypes2==0),cellrates(celltypes2==1),'method','approximate');
    fprintf('Unresponsive vs responsive z= %.2f, p=%.4f \n',stats.zval,p);
    [p,h,stats]=ranksum(cellrates(celltypes2==0),cellrates(celltypes2==2),'method','approximate');
    fprintf('Unresponsive vs selective z= %.2f, p=%.4f \n',stats.zval,p);
    [p,h,stats]=ranksum(cellrates(celltypes2==1),cellrates(celltypes2==2),'method','approximate');
    fprintf('responsive vs Selective z= %.2f, p=%.4f \n',stats.zval,p);
    fprintf('\n');
    % sgtitle doesnt really make sense here, gotta add numbers in xlabels
    sgtitle(sprintf('Place Fields of %s cells from %s',type, region{i}));
end
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

figure; scatter(objscores(celltypes==0),trajcorrs(celltypes==0),12,'filled');
hold on;
scatter(objscores(celltypes>0),trajcorrs(celltypes>0),12,'filled');
xlabel('object responsiveness score (chg from baseline)'); ylabel('outbound pf correlation');
legend('obj nonresponsive','obj responsive');
title(sprintf('%s %s cells',region{i}, type)); box off;
% first filter based on has a trajectory, then run corr on trajectories


% Does odor selectivity match trajectory selectivity (mean rates across
% runs)
% this is signed dprime effect size
odordiffs=cellfun(@(a) a(4), {cellPool.OdorSelective}).*cellfun(@(a) (((a(1)-a(2))>0)-.5)*2, {cellPool.OdorMeans});
% this is dprime effect size betweeen peak rates across runs
trajdiffs=cellfun(@(a) dprime(a{1},a{2}).*((mean(a{1})>mean(a{2}))-.5)*2, {cellPool.RunRates});

figure; scatter(odordiffs(celltypes<2),trajdiffs(celltypes<2),12,'filled');
hold on;
scatter(odordiffs(celltypes==2),trajdiffs(celltypes==2),12,'filled');
xlabel('d'' obj selectivity'); ylabel('d'' outbout pf peak rate');
legend('obj nonselective','obj selective');
title(sprintf('%s %s cells',region{i}, type));




% Does odor selectivity match diff in peak rates fo each run
% signed d' of odor coding
odordiffs=cellfun(@(a) a(4), {cellPool.OdorSelective}).*cellfun(@(a) (((a(1)-a(2))>0)-.5)*2, {cellPool.OdorMeans});
% this is just a raw difference in max
trajdiffs=cellfun(@(a) (max(a{1})-max(a{2}))/...
    (max(a{1})+max(a{2})), {cellPool.LinPlaceFields});


figure; s=scatter(odordiffs(celltypes<2),trajdiffs(celltypes<2),12,'filled');
hold on;
s(2)=scatter(odordiffs(celltypes==2),trajdiffs(celltypes==2),12,'filled');
xlabel('Odor Selectivity (d'')'); ylabel('Outbound Trajectory Selectivity (grand mean diff/sum)');
title(sprintf('%s %s cells',region{i}, type));
hold on;
plot([-1 1],[0 0],'k'); plot([0 0],[-1 1],'k');
legend(s,'Non Selective','Odor Selective'); box off;
xlim([-1 1]); ylim([-1 1]);
[rho,p]=corr(odordiffs(celltypes<2)',trajdiffs(celltypes<2)','rows','complete');
[rho2,p2]=corr(odordiffs(celltypes==2)',trajdiffs(celltypes==2)','rows','complete');
fprintf('for full %s pop r^2=%.2f, p=%.4f, for selective cells, r^2=%.2f, p=%.4f \n',...
    region{i},rho,p,rho2,p2);


% howabout whether the cell had a place field ont he same or opposing
% trajectory


%somethiing like none, same, opposite both
% we need object selectivity, and the direction
% thats objSel,  now we need pf presence, so it will have to be something
% like 0, 1, or 2, and then if it has one, we multiply it by the diff
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
legend('Object Responsive','Object Selective');
ylabel('% of units');
xtips = b(1).XEndPoints; ytips = b(1).YEndPoints;
labels = {sum(pfnums(objSel~=0)==0),...
     sum(pfnums(objSel~=0)==objSel(objSel~=0)),...
     sum(pfnums(objSel~=0)==-objSel(objSel~=0)),...
     sum(pfnums(objSel~=0)==2)};
text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');

xtickangle(340); box off;
end


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

for i=1:length(region)
% filter by region and celltype and has to fire during at least one run
cellfilt=cellfun(@(a) contains(a,region{i},'IgnoreCase',true),{SuperUnits.area}) &...
    cellfun(@(a) contains(a,type,'IgnoreCase',true),{SuperUnits.type}) &...
    cellfun(@(a) any(a),{SuperUnits.FiresDuringRun});

% of the splitter cells, whats their object selectivity?
cellPool=SuperUnits(cellfilt);
Spllitters=cellfun(@(a) a(1)>0 & a(2)<.05, {cellPool.SplitterScore});

% now compare rates of odor responsiveness and odor selectivity of the
% splitters to the non splitters



end


