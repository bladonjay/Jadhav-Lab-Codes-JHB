%% PlotObjPlaceCells



%% this plots all the object selective cells
% it does so for each day individually, not for each epoch in a day

%%
Params.region='PFC';
Params.nameappend='-PlaceSelective';
Params.speedthresh=3;
Params.timesmooth=8; % bins
Params.maxtimejump=1;
saveoutdir=uigetdir('e:\',sprintf('Region is %s, name to append is %s',Params.region,Params.nameappend));
%%

for ses=1:length(SuperRat)
    
    %----figure out which cells we want to use
    % cells in the region that code the right thing
    %CodingCell=cellfun(@(a) a(4)<.05 & a(4)>0, {SuperRat(ses).units.OdorResponsive});
    %CodingCell=cellfun(@(a) a(3)==1, {SuperRat(ses).units.OdorSelective});
    CodingCell=cellfun(@(a) any(a==1), {SuperRat(ses).units.PFexist});                      %***** CHANGE HERE

    % NO INS
    inRegion=cellfun(@(a) contains(a,Params.region), {SuperRat(ses).units.area});
    isPyram=cellfun(@(a) contains(a,'pyr'),{SuperRat(ses).units.type}); % already classified
    UseCells=find(CodingCell & inRegion & isPyram);
    
    %---- get trial data
    % get curves for this sessions trials (zdefunct, data are
    % already in struct
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    fulltrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10,...
        trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
    % correct trials, and ID their odor
    oktrials=fulltrialmat(:,4)==1 & ismember(fulltrialmat(:,5),SuperRat(ses).RunEpochs);
    trialmat=fulltrialmat(oktrials,:); % no learning data!!!

    %---- get tracking lined up (remove blocks with no spikes tho...)
    %  use only spikes and tracking where the animal is moving
    % were doign this differently because the lincoords is already
    % stitched (so there are major time jumps)
    % have to smooth velocity data because it oscillates everty 3 time
    % indices
    tracking=SuperRat(ses).tracking.data;
    % pull only run data
    okepoch=ismember(tracking(:,6),SuperRat(ses).RunEpochs);
    disptracking=tracking(okepoch,:); % this is all the run data
    % 8 bins to smooth velocity, speed less than 3 ditch
    tooslow=SmoothMat2(disptracking(:,5),[0 50],Params.timesmooth)<=Params.speedthresh;
    
    % this scrubs the slow data to calculate occupancy
    thistrack=disptracking(~tooslow,:);
    % this removes spikes inbetween run epochs
    jumps=find(diff(thistrack(:,1))>Params.maxtimejump);
    goodepochs=[[thistrack(1); thistrack(jumps+1,1)] [thistrack(jumps,1); thistrack(end,1)]];

    % plot each cell
    for i=1:length(UseCells)
        figure;
        spikedata=SuperRat(ses).units(UseCells(i)).ts(:,1);
        
        % remove spikes between the tracking
        subplot(4,1,1);
        thesespikes=EpochCoords(spikedata,goodepochs);
        fprintf('%d of %d spikes are during runs \n',length(thesespikes),length(spikedata));
        plot(disptracking(:,2),disptracking(:,3),'k');
        hold on;
        spikex=interp1(thistrack(:,1),thistrack(:,2),thesespikes);
        spikey=interp1(thistrack(:,1),thistrack(:,3),thesespikes);
        plot(spikex,spikey,'r.'); box off;
        
        title([SuperRat(ses).name ' ' num2str(ses) ' ' SuperRat(ses).units(UseCells(i)).area,...
            ' ' SuperRat(ses).units(UseCells(i)).type ' ' num2str(UseCells(i))]);
        %end
        axis tight; box off; axis off;
        
        % now make place plot here
        subplot(4,1,2);
        session=struct('edit_coords',thistrack);
        thisunit=struct('ts',thesespikes);
        % this is the really basic open field place plotting script
        [ratemap,~,finalcolormap]=cell_SmoothPlacePlot(session,thisunit,...
            'Factor',2,'suppress',1,'gaussdev',1.5,'ColorScheme','parula');
        image(finalcolormap); set(gca,'YDir','normal');
        title(sprintf('         Max: %.2f Hz',max(linearize(ratemap))));
        box off; axis off;

        % going to do a box and whisker plot for odor rates
        subplot(4,1,3);
        
        [~,spkevs,~,trspkinds,~,evspikes]=event_spikes(spikedata,trialmat(:,1),...
            0,trialmat(:,2)-trialmat(:,1));
        spknums=cellfun(@(a) length(a), trspkinds);
        % scrub any block that is zero
        spikesperblock=accumarray(trialmat(:,5),spknums);
        trialsperblock=accumarray(trialmat(:,5),1);
        %find those inds and remove them
        keeptrials=sum(find(spikesperblock>trialsperblock/2)'==trialmat(:,5),2)>0;
        odorid=trialmat(:,3);

        %or you can do the raster for time locked firing rate
        subplot(4,1,3);
          [~,~,~,~,~,evspikes]=event_spikes(spikedata,trialmat(keeptrials,1),...
            2,trialmat(keeptrials,2)-trialmat(keeptrials,1)+2);
        thismat=trialmat(keeptrials,:);
         for q=1:length(evspikes)
             if ~isempty(evspikes{q})
                plot(linearize(repmat(evspikes{q}',3,1)),linearize(repmat([0 1 nan]',1,length(evspikes{q}))+q),'k');
                hold on;
                plot([thismat(q,2)-thismat(q,1) thismat(q,2)-thismat(q,1)],[q-1 q],'r');
             end
             ylim([0 length(evspikes)+1]);
         end

        % remember LR10
        
        if any(spkevs(keeptrials))
            
            boxScatterplot(spkevs(keeptrials),double(odorid(keeptrials)==0),'yLabel','Rate',...
                'xLabels',{sprintf('Odor Left %.2f Hz',SuperRat(ses).units(UseCells(i)).OdorMeans(1,1)),...
                sprintf('Odor Right %.2f Hz',SuperRat(ses).units(UseCells(i)).OdorMeans(1,2))},...
                'position',{},'plotBox',false);
            % now line bar
            odorl=spkevs(keeptrials & odorid==1); odorr=spkevs(keeptrials & odorid==0);
            hold on; plot(0,nanmean(odorl),'k.','MarkerSize',26);
            plot([0 0],[max([0 nanmean(odorl)-2*nanstd(odorl)]) nanmean(odorl)+2*nanstd(odorl)],'k.-','LineWidth',2);
            plot(1,nanmean(odorr),'k.','MarkerSize',26); hold on;
            plot([1 1],[max([0 nanmean(odorr)-2*nanstd(odorr)]) nanmean(odorr)+2*nanstd(odorr)],'k.-','LineWidth',2);
            set(gca,'XTick',[0 1],'XTickLabel',{sprintf('Odor Left %.2f Hz',SuperRat(ses).units(UseCells(i)).OdorMeans(1,1)),...
                sprintf('Odor Right %.2f Hz',SuperRat(ses).units(UseCells(i)).OdorMeans(1,2))});
            xlim([-0.5 1.5]); ylabel('Firing rate, Hz');
            title(sprintf('A=%.2f B=%.2f, p=%.3f', SuperRat(ses).units(UseCells(i)).OdorMeans(1,1),...
                SuperRat(ses).units(UseCells(i)).OdorMeans(1,2),...
                SuperRat(ses).units(UseCells(i)).OdorSelective(1,2)));
            box off;
        else
            title('No Spikes During odor delivery');
        end
        % now to make the linearized place plots
        % for runs its out left, out right in left in right (1:4)
        subplot(4,1,4);
        mycolors=lines(2);
        % out left
        if iscell(SuperRat(ses).units(UseCells(i)).LinPlaceFields)
            try
                p=plot(SuperRat(ses).units(UseCells(i)).LinPlaceFields{1}(1,:),'color',mycolors(1,:),'LineWidth',3); hold on;
            end
            try
                plot(-SuperRat(ses).units(UseCells(i)).LinPlaceFields{1}(2,:),'color',mycolors(1,:),'LineWidth',3);
            end
            % in left
            try
             p(2)=plot(SuperRat(ses).units(UseCells(i)).LinPlaceFields{1}(3,:),'color',mycolors(2,:),'LineWidth',3);
                % in right
                plot(-SuperRat(ses).units(UseCells(i)).LinPlaceFields{1}(4,:),'color',mycolors(2,:),'LineWidth',3);
            end
            plot([0 100],[0 0],'k'); % get a zero line in there
            legend(p,'Outbound','Inbound');
            oldlim=max(abs(get(gca,'Ylim'))); ylim([-oldlim oldlim]);
            box off; ylabel({'Rate, Hz';'Right Arm      Left Arm'});
        end
        
        
        %legend('Left Outbound','Right Outbound','Left Return','Right Return');
        pfstats=SuperRat(ses).units(UseCells(i)).FieldProps;
        [info,bestfield]=max(pfstats.info);
        title(sprintf('Max=%.1fHz, Info=%.2f, Width=%2.f',max(pfstats.PFmax),...
            max(pfstats.info),pfstats.PFsize(bestfield)));
        set(gca,'XTick',[0 50 100],'XTickLabel',{'Odor Port','Decision Point','Reward Port'});
        xlabel('Lineaized Position on maze');
        sgtitle(sprintf('Ses %s Unit %s # %d', SuperRat(ses).name,...
            SuperRat(ses).units(UseCells(i)).area,UseCells(i)));
        % set the plot size
        set(gcf,'Position',[-450 100 300 900]);
        % save out
        savefig(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(UseCells(i)).area,...
            '-' SuperRat(ses).units(UseCells(i)).type num2str(UseCells(i)) Params.nameappend]));
        saveas(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(UseCells(i)).area,...
            '-' SuperRat(ses).units(UseCells(i)).type num2str(UseCells(i)) Params.nameappend]),'tif');
        close(gcf);
    end
end
%%


% get a superrat of units
SuperUnits=orderfields(SuperRat(1).units);
SuperUnits=rmfield(SuperUnits,'csi');
for i=2:length(SuperRat)
    theseunits=orderfields(SuperRat(i).units);
    % for some reason csi is in some of these structs
    if isfield(theseunits,'csi'), theseunits=rmfield(theseunits,'csi'); end
    SuperUnits=[SuperUnits theseunits];
end

%%

location='CA1';
%location='PFC';

% whose lin place field plots have information?
% lets just aggregate all the odor selective cells and get the % of them
% who have spatial information, and the % of them who have place fields

% i think it has to have some minimum # of spikes... maybe has to fire on
% at least one run or has to fire during odor period?
ActivePop=cellfun(@(a) any(a), {SuperUnits.FiresDuringRun}) |...
    cellfun(@(a) a(1)>0, {SuperUnits.OdorResponsive});

ispyram=cellfun(@(a) contains(a,'pyr'), {SuperUnits.type});


%ActivePop=cell2mat(cellfun(@(a) any(a(:)), {SuperUnits.FiresDuringRun},'UniformOutput', false));
isHere=cell2mat(cellfun(@(a) contains(a,location), {SuperUnits.area},'UniformOutput', false));

% it may be good to try to use a table


usecell=ActivePop(:) & isHere(:) & ispyram(:);

% 1. which are obj selective cells
allstats=cell2mat(cellfun(@(a) a(:,3), {SuperUnits(usecell).OdorSelective}, 'UniformOutput', false)');
% 2. which have place fields? (out and return)
allstats(:,2)=cell2mat(cellfun(@(a) sum(a(1:4)), {SuperUnits(usecell).PFexist}, 'UniformOutput', false)');
% 3. Mean rate during runs
allstats(:,3)=cell2mat(cellfun(@(a) mean(a{1}(:)), {SuperUnits(usecell).LinPlaceFields}, 'UniformOutput', false)');
% 4. information in bits
allstats(:,4)=cell2mat(cellfun(@(a) max(a.info,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 5. sparsity
allstats(:,5)=cell2mat(cellfun(@(a) max(a.sparsity,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 6. place field size
allstats(:,6)=cell2mat(cellfun(@(a) min(a.PFsize,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 7. peak place field rate
allstats(:,7)=cell2mat(cellfun(@(a) max(a.PFmax,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 8. Trajectory score (diff over sum for runs)
allstats(:,8)=cell2mat(cellfun(@(a) a(:,1), {SuperUnits(usecell).TrajScores}, 'UniformOutput', false)');
% 9 odor selectivity score
allstats(:,9)=cell2mat(cellfun(@(a) a(:,1), {SuperUnits(usecell).OdorSelective}, 'UniformOutput', false)');
% 10 which object is it selective to
allstats(:,10)=cell2mat(cellfun(@(a) a(1)>0, {SuperUnits(usecell).OdorSelective}, 'UniformOutput', false)');
% 11 is it trajectory selective (pvalue)
allstats(:,11)=cell2mat(cellfun(@(a) a(:,2), {SuperUnits(usecell).TrajScores}, 'UniformOutput', false)');
% 12 difference between peaks as % of sum... but only of they're above say 2 Hz
allstats(:,12)=cellfun(@(a) (max(a{1}(1,:))-max(a{1}(2,:)))/(max(a{1}(1,:))+max(a{1}(2,:))),...
    {SuperUnits(usecell).LinPlaceFields});
% 13 is the index of that difference
allstats(:,13)=cell2mat(cellfun(@(a) find(a{1}(1,:)==max(a{1}(1,:)),1,'first')-find(a{1}(2,:)==max(a{1}(2,:)),1,'first'),...
    {SuperUnits(usecell).LinPlaceFields}, 'UniformOutput', false));
% 14 is it odor responsive? (elevated rate)
allstats(:,14)=cellfun(@(a) a(2)>0 && a(4)<.05, {SuperUnits(usecell).OdorResponsive});
% 15 is the correlation between trajectories 1 and 2
allstats(:,15)=cell2mat(cellfun(@(a) corr(a{1}(1,:)',a{1}(2,:)'), {SuperUnits(usecell).LinPlaceFields}, 'UniformOutput', false))';
% 16 is does it fire on each trajectory? (already has to fire +2 on a
% trajectory
allstats(:,16)=cell2mat(cellfun(@(a) sum(a(1:2)), {SuperUnits(usecell).FiresDuringRun}, 'UniformOutput', false)');
% 17 whats the change in rate from before to odor sampling (zscored)
% (a-b)/Std(a) % so kindof
allstats(:,17)=cellfun(@(a) a(2)/a(3), {SuperUnits(usecell).OdorResponsive});
% 18 has a peak fr of above .5 hz on any trajectory
allstats(:,18)=cellfun(@(a) sum(a.PFmax>0.5), {SuperUnits(usecell).FieldProps});

allstats(:,19)=cell2mat({SuperUnits(usecell).meanrate});

allstats(:,20)=cellfun(@(a) a(1), {SuperUnits(usecell).OdorResponsive});

%%
% does the rate during running relate to the rate during odor sampling
% it seems that the robustly firing neurons during run, depress their rate
% during odor sampling, and the opposite is true as well 
figure; subplot(2,3,1);
scatter(log10(allstats(:,7)),allstats(:,17),5,'filled');
xlabel('Peak run rate(Log_1_0)'); ylabel('Object responsivity (Norm \Delta rate)');
okcells=~isnan(allstats(:,7)) & ~isnan(allstats(:,17));
[a,b]=corr(allstats(okcells,7),allstats(okcells,17),'Type','Pearson');
title(sprintf('R^2= %.2f, P=%.4f',a,b));

subplot(2,3,2);
% is equivalnce across runs related to odor rate?

scatter(allstats(:,17),allstats(:,15),5,'filled');
xlabel('Object resposiveness ( \Delta rate)'); ylabel('Correlation between trajectories');
okcells=~isnan(allstats(:,17)) & ~isnan(allstats(:,15));
[a,b]=corr(allstats(okcells,17),allstats(okcells,15),'Type','Pearson');
title(sprintf('R^2= %.2f, P=%.4f',a,b));

% is trajectory equivalence just a factor of rate?
subplot(2,3,3);
scatter(log10(allstats(:,7)),allstats(:,15),5,'filled');
ylabel('Correlation between trajectories'); xlabel('Peak Run Rate (Log_1_0) Hz');
okcells=~isnan(allstats(:,7)) & ~isnan(allstats(:,15));
[a,b]=corr(allstats(okcells,7),allstats(okcells,15),'Type','Pearson');
title(sprintf('R^2= %.2f, P=%.4f',a,b));

% okay does this just dependn on overall rate????
subplot(2,3,4);
%okcells=~isnan(allstats(:,3)) & ~isnan(allstats(:,19)) & ~isnan(allstats(:,20));
%scatter(log10(allstats(okcells,3)),log10(allstats(okcells,20)),25,log10(allstats(okcells,19)),'filled');
%xlabel('Run Rate'); ylabel('Odor Rate');

okcells=~isnan(allstats(:,3)) & ~isnan(allstats(:,9)) & ~isnan(allstats(:,14));
scatter(log10(allstats(okcells,3)),log10(allstats(okcells,20)),25,log10(allstats(okcells,9)),'filled');
xlabel('Mean Run Rate'); ylabel('Mean Odor Rate');
legend('Color is odor selectivity (diff/sum)')

subplot(2,3,5);


okcells=~isnan(allstats(:,7)) & ~isnan(allstats(:,17)) & ~isnan(allstats(:,14));
scatter(log10(allstats(okcells,7)),allstats(okcells,17),25,log10(allstats(okcells,9)),'filled');
xlabel('Peak PF rate'); ylabel('Odor Responsivity (Norm \Delta Rate)');
legend('Color is odor selectivity (diff/sum)')

subplot(2,3,6);
okcells=~isnan(allstats(:,7)) & ~isnan(allstats(:,8)) & ~isnan(allstats(:,14));
scatter(log10(allstats(okcells,7)),allstats(okcells,8),25,allstats(okcells,1)/5,'filled');
xlabel('Peak PF rate'); ylabel('Trajectory selectivity (Norm \Delta Rate)');
legend('Color is odor selectivity (yellow is selective)')

% what proportion of object selective cells fire on the maze, what
% proportion of object responsive cells fire on the maze (have a peak >3 Hz
% how about one field vs more than one

barx=[nanmean(allstats(allstats(:,1)==1,18)==1), nanmean(allstats(allstats(:,1)==1,18)>1);...
    nanmean(allstats(allstats(:,14)==1,18)==1), nanmean(allstats(allstats(:,14)==1,18)>1)];
figure; sp=subplot(1,2,1);
bar(barx,'stacked');
set(gca,'XTickLabel',{'Object Selective','Object Responsive'});
ylabel('Proportion of cells with a track peak >0.5 Hz');
legend('Fires on one trajectory','more than one trajectory');

% howabout has pf

barx=[nanmean(allstats(allstats(:,1)==1,2)==1), nanmean(allstats(allstats(:,1)==1,2)>1);...
    nanmean(allstats(allstats(:,14)==1,2)==1), nanmean(allstats(allstats(:,14)==1,2)>1)];
sp(2)=subplot(1,2,2);
bar(barx,'stacked');
set(gca,'XTickLabel',{'Object Selective','Object Responsive'});
ylabel('Proportion of cells with a track field');
legend('field on one trajectory','field on multiple trajectories');
linkaxes(sp,'y');
%
% okay now lets look at the place field characteristics of these cells
% first I want to check that the pf peaks and sizes are the same
% but if the cell has multiple place fields, which do i take?
% first lets take all of them


% start with place field characteristics and assign cell afterwards
allpfs=[];

for i=1:length(SuperUnits)% do th
    if usecell(i) % only hc non ins
        % pfmax, pfsze, pf sparsity
        thispf=[SuperUnits(i).FieldProps.PFmax' SuperUnits(i).FieldProps.PFsize' SuperUnits(i).FieldProps.sparsity'];
        thispf(:,7)=1:4;
        % has to have a place field to analyze its spatial properties
        thispf=thispf(SuperUnits(i).PFexist,:);
        % is it odor selective?
        thispf(:,4)=SuperUnits(i).OdorSelective(3)==1;
        % sig elevated rate
        thispf(:,5)=SuperUnits(i).OdorResponsive(4)<.05 & SuperUnits(i).OdorResponsive(2)>0;
        % what is the odor preference?
        thispf(:,6)=SuperUnits(i).OdorSelective(1);
        % which route is this?
        
        allpfs=[allpfs;thispf];
    end
end
% allpfs is 1. pfmax, 2 pfsize 3 sparsity, 4 is cell odor selective 5 is
% cell odor responsive 6 is direction of odor selectivity
% nan out bad cells

figure; 
subplot(3,2,1);
boxScatterplot(allpfs(:,1),allpfs(:,4),'position',[],'Xlabels',{'not','Odor Selective'});
ylabel('PF peak rate'); 
title(sprintf('nonselective: %.2f selective %.2f p=%.4f',nanmean(allpfs(allpfs(:,4)==0,1)),...
    nanmean(allpfs(allpfs(:,4)==1,1)),ranksum(allpfs(allpfs(:,4)==0,1),allpfs(allpfs(:,4)==1,1))));

subplot(3,2,2);
boxScatterplot(allpfs(:,1),allpfs(:,5),'position',[],'Xlabels',{'not','Odor Responsive'});
ylabel('PF peak rate');
title(sprintf('nonresponsive: %.2f responsive %.2f p=%.4f',nanmean(allpfs(allpfs(:,5)==0,1)),...
    nanmean(allpfs(allpfs(:,5)==1,1)),ranksum(allpfs(allpfs(:,5)==0,1),allpfs(allpfs(:,5)==1,1))));

subplot(3,2,3);
boxScatterplot(allpfs(:,2),allpfs(:,4),'position',[],'Xlabels',{'not','Odor Selective'});
ylabel('PF size'); 
title(sprintf('nonselective: %.2f selective %.2f p=%.4f',nanmean(allpfs(allpfs(:,4)==0,2)),...
    nanmean(allpfs(allpfs(:,4)==1,2)),ranksum(allpfs(allpfs(:,4)==0,2),allpfs(allpfs(:,4)==1,2))));

subplot(3,2,4);
boxScatterplot(allpfs(:,2),allpfs(:,5),'position',[],'Xlabels',{'not','Odor Responsive'});
ylabel('PF size');
title(sprintf('nonresponsive: %.2f responsive %.2f p=%.4f',nanmean(allpfs(allpfs(:,5)==0,2)),...
    nanmean(allpfs(allpfs(:,5)==1,2)),ranksum(allpfs(allpfs(:,5)==0,2),allpfs(allpfs(:,5)==1,2))));

subplot(3,2,5);
boxScatterplot(log(allpfs(:,3)),allpfs(:,4),'position',[],'Xlabels',{'not','Odor Selective'});
ylabel('sparsity'); 
title(sprintf('nonselective: %.2f selective %.2f p=%.4f',nanmean(allpfs(allpfs(:,4)==0,3)),...
    nanmean(allpfs(allpfs(:,4)==1,3)),ranksum(allpfs(allpfs(:,4)==0,3),allpfs(allpfs(:,4)==1,3))));

subplot(3,2,6);
boxScatterplot(log(allpfs(:,3)),allpfs(:,5),'position',[],'Xlabels',{'not','Odor Responsive'});
ylabel('PF sparsity');
title(sprintf('nonresponsive: %.2f responsive %.2f p=%.4f',nanmean(allpfs(allpfs(:,5)==0,3)),...
    nanmean(allpfs(allpfs(:,5)==1,3)),ranksum(allpfs(allpfs(:,5)==0,3),allpfs(allpfs(:,5)==1,3))));



%% so lets zoom in here, do odor cells show the same discrimination across runs?

% basically, you would expect the odor ensemble to not interfere with the
% run ensemble??

% so i guess the question is whether a odor cell is apt to fire on the
% opposite run?


%% 1. do cells that respond to odors respond to place differently



