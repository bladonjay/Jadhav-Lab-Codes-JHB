%% PlotObjPlaceCells_old


% for ses=1:length(SuperRat)

% what proportion of cells are odor selective
cumOdor={};
for i=1:length(SuperRat)
    cumobj{i}=cell2mat({SuperRat(i).units.OdorSelective}');
end
temp=cell2mat(cumobj');
fprintf(' %d cells, or %.2f %% of all cells coded odors \n',...
    nansum(temp(:,3)), nanmean(temp(:,3)*100));

%%
saveoutdir=uigetdir;
nameappend='-objplaceplot';

%% prelim plot for all the odor selective cells (if we're plotting each run)
for ses=1:length(SuperRat)
    % odor selective cells
    objectcells=find(cellfun(@(a) any(a(:,3)==1), {SuperRat(ses).units.OdorSelective}));
    trialdata=SuperRat(ses).trialdata;
    % plot each cell
    for i=1:length(objectcells)
        figure;
        Runs=SuperRat(ses).RunEpochs;
        tracking=SuperRat(ses).tracking;
        spikedata=SuperRat(ses).units(objectcells(i)).ts;
        trialmat=[trialdata.sniffstart trialdata.sniffend...
            trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
        
        % for each day, do a top stick and dot plot, and a bottom place
        % plot
        for k=1:length(Runs)
            % top two are place maps
            thistrack=tracking.data(tracking.data(:,6)==Runs(k),:);
            thesespikes=spikedata(spikedata>thistrack(1,1) & spikedata<thistrack(end,1));
            subplot(4,length(Runs),k);
            plot(thistrack(:,2),thistrack(:,3),'k');
            hold on;
            spikex=interp1(thistrack(:,1),thistrack(:,2),thesespikes);
            spikey=interp1(thistrack(:,1),thistrack(:,3),thesespikes);
            plot(spikex,spikey,'r*'); box off;
            %legend('trajectory','Spikes');
            if k==1
                title([SuperRat(ses).name ' ' num2str(ses) ' ' SuperRat(ses).units(objectcells(i)).area,...
                    ' ' SuperRat(ses).units(objectcells(i)).type ' ' num2str(objectcells(i))]);
            end
            axis tight; box off; axis off;
            
            % now make place plot here
            subplot(4,length(Runs),k+length(Runs));
            session=struct('edit_coords',tracking.data(tracking.data(:,6)==Runs(k),1:3));
            
            % this is the really basic open field place plotting script
            [ratemap,~,finalcolormap]=cell_SmoothPlacePlot(session,SuperRat(ses).units(objectcells(i)),...
                'Factor',3,'suppress',1,'gaussdev',1.5);
            image(finalcolormap); set(gca,'YDir','normal');
            title(sprintf('Place Plot Max: %.2f',max(linearize(ratemap))));
            box off; axis off;
            
            % get curves for this sessions trials (defunct, data are
            % already in struct
            goodtrials=trialmat(:,4)==1 & trialmat(:,5)==Runs(k);
            odorid=trialmat(goodtrials,3);
            %[~,~,~,~,~,spikets]=event_spikes(thesespikes,trialmat(goodtrials,1),...
            %    abs(timeedges(1)),timeedges(end));
            
            % going to do a box and whisker plot for odor rates
            [~,spkevs]=event_spikes(thesespikes,trialmat(goodtrials,1),...
                0,trialmat(goodtrials,2)-trialmat(goodtrials,1));
            % remember LR10
            subplot(4,length(Runs),k+length(Runs)*2);
            boxScatterplot(spkevs,double(trialmat(goodtrials,3)==0),'yLabel','Rate',...
                'xLabels',{'Odor left','Odor Right'},'position',{});
            title(sprintf('A=%.2f B=%.2f, p=%.3f', SuperRat(ses).units(objectcells(i)).OdorRates(k,1),...
                SuperRat(ses).units(objectcells(i)).OdorRates(k,2),...
                SuperRat(ses).units(objectcells(i)).OdorSelective(k,2)));
            box off;
            
            % now to make the linearized place plots
            % for runs its out left, out right in left in right (1:4)
            subplot(4,length(Runs),k+length(Runs)*3);
            plot(SuperRat(ses).units(objectcells(i)).LinPlaceFields{k}(1,:)); hold on;
            plot(-SuperRat(ses).units(objectcells(i)).LinPlaceFields{k}(2,:));
            plot(SuperRat(ses).units(objectcells(i)).LinPlaceFields{k}(3,:));
            plot(-SuperRat(ses).units(objectcells(i)).LinPlaceFields{k}(4,:));
            plot([0 100],[0 0],'k'); % get a zero line in there
            
            box off; ylabel('Rate, Hz');
            %legend('Left Outbound','Right Outbound','Left Return','Right Return');
            pfstats=SuperRat(ses).units(objectcells(i)).FieldProps;
            title(sprintf('Max=%.1fHz, Info=%.2f, Width=%2.f',max(pfstats.PFmax(k,:)),...
                max(pfstats.info(k,:)),min(pfstats.PFsize(k,:))));
            set(gca,'XTick',[0 100],'XTickLabel',{'Odor Port','Reward Port'});
        end
        % name the unit
        sgtitle(sprintf('Ses %s Unit %s %s # %d', SuperRat(ses).name,...
            SuperRat(ses).units(objectcells(i)).area,...
            SuperRat(ses).units(objectcells(i)).type,objectcells(i)));
        % set the plot size
        set(gcf,'Position',[250 250 300*length(Runs) 800]);
        % save out
        savefig(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(objectcells(i)).area,...
            '-' SuperRat(ses).units(objectcells(i)).type num2str(objectcells(i)) nameappend]));
        saveas(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(objectcells(i)).area,...
            '-' SuperRat(ses).units(objectcells(i)).type num2str(objectcells(i)) nameappend]),'tif');
        close(gcf);
    end
end
%% this plots all the object selective cells
% it does so for each day individually, not for each epoch in a day

saveoutdir=uigetdir;
nameappend='-objplaceplotObjCell';

for ses=1:length(SuperRat)
    % odor selective cells
    objectcells=find(cellfun(@(a) any(a(:,3)==1), {SuperRat(ses).units.OdorSelective}));
    trialdata=SuperRat(ses).trialdata;
    % plot each cell
    for i=1:length(objectcells)
        figure;
        tracking=SuperRat(ses).tracking;
        thistrack=tracking.data;
        spikedata=SuperRat(ses).units(objectcells(i)).ts;
        trialmat=[trialdata.sniffstart trialdata.sniffend...
            trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
        % only take useful spikes
        thesespikes=spikedata(spikedata>thistrack(1,1) & spikedata<thistrack(end,1));
        subplot(4,1,1);
        plot(thistrack(:,2),thistrack(:,3),'k');
        hold on;
        spikex=interp1(thistrack(:,1),thistrack(:,2),thesespikes);
        spikey=interp1(thistrack(:,1),thistrack(:,3),thesespikes);
        plot(spikex,spikey,'r*'); box off;
        %legend('trajectory','Spikes');
        title([SuperRat(ses).name ' ' num2str(ses) ' ' SuperRat(ses).units(objectcells(i)).area,...
            ' ' SuperRat(ses).units(objectcells(i)).type ' ' num2str(objectcells(i))]);
        axis tight; box off; axis off;
        
        % now make place plot here
        subplot(4,1,2);
        session=struct('edit_coords',tracking.data(:,1:3));
        
        % this is the really basic open field place plotting script
        [ratemap,~,finalcolormap]=cell_SmoothPlacePlot(session,SuperRat(ses).units(objectcells(i)),...
            'Factor',3,'suppress',1,'gaussdev',1.5);
        image(finalcolormap); set(gca,'YDir','normal');
        title(sprintf('Place Plot Max: %.2f',max(linearize(ratemap))));
        box off; axis off;
        
        % get curves for this sessions trials (defunct, data are
        % already in struct
        goodtrials=trialmat(:,4)==1;
        odorid=trialmat(goodtrials,3);
        %[~,~,~,~,~,spikets]=event_spikes(thesespikes,trialmat(goodtrials,1),...
        %    abs(timeedges(1)),timeedges(end));
        
        % going to do a box and whisker plot for odor rates
        [~,spkevs]=event_spikes(thesespikes,trialmat(goodtrials,1),...
            0,trialmat(goodtrials,2)-trialmat(goodtrials,1));
        % remember LR10
        subplot(4,1,3);
        boxScatterplot(spkevs,double(trialmat(goodtrials,3)==0),'yLabel','Rate',...
            'xLabels',{'Odor left','Odor Right'},'position',{});
        title(sprintf('A=%.2f B=%.2f, p=%.3f', SuperRat(ses).units(objectcells(i)).OdorMeans(1,1),...
            SuperRat(ses).units(objectcells(i)).OdorMeans(1,2),...
            SuperRat(ses).units(objectcells(i)).OdorSelective(1,2)));
        box off;
        
        % now to make the linearized place plots
        % for runs its out left, out right in left in right (1:4)
        subplot(4,1,4);
        plot(SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(1,:)); hold on;
        plot(-SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(2,:));
        plot(SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(3,:));
        plot(-SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(4,:));
        plot([0 100],[0 0],'k'); % get a zero line in there
        
        box off; ylabel('Rate, Hz');
        %legend('Left Outbound','Right Outbound','Left Return','Right Return');
        pfstats=SuperRat(ses).units(objectcells(i)).FieldProps;
        title(sprintf('Max=%.1fHz, Info=%.2f, Width=%2.f',max(pfstats.PFmax),...
            max(pfstats.info),min(pfstats.PFsize)));
        set(gca,'XTick',[0 100],'XTickLabel',{'Odor Port','Reward Port'});
        sgtitle(sprintf('Ses %s Unit %s %s # %d', SuperRat(ses).name,...
            SuperRat(ses).units(objectcells(i)).area,...
            SuperRat(ses).units(objectcells(i)).type,objectcells(i)));
        % set the plot size
        set(gcf,'Position',[250 250 300*length(Runs) 800]);
        % save out
        savefig(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(objectcells(i)).area,...
            '-' SuperRat(ses).units(objectcells(i)).type num2str(objectcells(i)) nameappend]));
        saveas(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(objectcells(i)).area,...
            '-' SuperRat(ses).units(objectcells(i)).type num2str(objectcells(i)) nameappend]),'tif');
        close(gcf);
    end
end
%% This plots all the cells with   'place fields'

saveoutdir=uigetdir;
nameappend='-objplaceplotPlaceCell';

for ses=1:length(SuperRat)
    % odor selective cells
    objectcells=find(cellfun(@(a) any(a(:,1:2)==1), {SuperRat(ses).units.PFexist}));
    trialdata=SuperRat(ses).trialdata;
    % plot each cell
    for i=1:length(objectcells)
        figure;
        tracking=SuperRat(ses).tracking;
        thistrack=tracking.data;
        spikedata=SuperRat(ses).units(objectcells(i)).ts;
        trialmat=[trialdata.sniffstart trialdata.sniffend...
            trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
        % only take useful spikes
        thesespikes=spikedata(spikedata>thistrack(1,1) & spikedata<thistrack(end,1));
        subplot(4,1,1);
        plot(thistrack(:,2),thistrack(:,3),'k');
        hold on;
        spikex=interp1(thistrack(:,1),thistrack(:,2),thesespikes);
        spikey=interp1(thistrack(:,1),thistrack(:,3),thesespikes);
        plot(spikex,spikey,'r*'); box off;
        %legend('trajectory','Spikes');
        %if k==1
        title([SuperRat(ses).name ' ' num2str(ses) ' ' SuperRat(ses).units(objectcells(i)).area,...
            ' ' SuperRat(ses).units(objectcells(i)).type ' ' num2str(objectcells(i))]);
        %end
        axis tight; box off; axis off;
        
        % now make place plot here
        subplot(4,1,2);
        session=struct('edit_coords',tracking.data(:,1:3));
        
        % this is the really basic open field place plotting script
        [ratemap,~,finalcolormap]=cell_SmoothPlacePlot(session,SuperRat(ses).units(objectcells(i)),...
            'Factor',3,'suppress',1,'gaussdev',1.5);
        image(finalcolormap); set(gca,'YDir','normal');
        title(sprintf('Place Plot Max: %.2f',max(linearize(ratemap))));
        box off; axis off;
        
        % get curves for this sessions trials (zdefunct, data are
        % already in struct
        goodtrials=trialmat(:,4)==1;
        odorid=trialmat(goodtrials,3);
        %[~,~,~,~,~,spikets]=event_spikes(thesespikes,trialmat(goodtrials,1),...
        %    abs(timeedges(1)),timeedges(end));
        
        % going to do a box and whisker plot for odor rates
        [~,spkevs]=event_spikes(thesespikes,trialmat(goodtrials,1),...
            0,trialmat(goodtrials,2)-trialmat(goodtrials,1));
        % remember LR10
        subplot(4,1,3);
        boxScatterplot(spkevs,double(trialmat(goodtrials,3)==0),'yLabel','Rate',...
            'xLabels',{'Odor left','Odor Right'},'position',{});
        title(sprintf('A=%.2f B=%.2f, p=%.3f', SuperRat(ses).units(objectcells(i)).OdorMeans(1,1),...
            SuperRat(ses).units(objectcells(i)).OdorMeans(1,2),...
            SuperRat(ses).units(objectcells(i)).OdorSelective(1,2)));
        box off;
        
        % now to make the linearized place plots
        % for runs its out left, out right in left in right (1:4)
        subplot(4,1,4);
        plot(SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(1,:)); hold on;
        plot(-SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(2,:));
        plot(SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(3,:));
        plot(-SuperRat(ses).units(objectcells(i)).LinPlaceFields{1}(4,:));
        plot([0 100],[0 0],'k'); % get a zero line in there
        
        box off; ylabel('Rate, Hz');
        %legend('Left Outbound','Right Outbound','Left Return','Right Return');
        pfstats=SuperRat(ses).units(objectcells(i)).FieldProps;
        title(sprintf('Max=%.1fHz, Info=%.2f, Width=%2.f',max(pfstats.PFmax),...
            max(pfstats.info),min(pfstats.PFsize)));
        set(gca,'XTick',[0 100],'XTickLabel',{'Odor Port','Reward Port'});
        sgtitle(sprintf('Ses %s Unit %s %s # %d', SuperRat(ses).name,...
            SuperRat(ses).units(objectcells(i)).area,...
            SuperRat(ses).units(objectcells(i)).type,objectcells(i)));
        % set the plot size
        set(gcf,'Position',[250 250 300 800]);
        % save out
        savefig(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(objectcells(i)).area,...
            '-' SuperRat(ses).units(objectcells(i)).type num2str(objectcells(i)) nameappend]));
        saveas(gcf,fullfile(saveoutdir,[SuperRat(ses).name '-' num2str(ses) '-' SuperRat(ses).units(objectcells(i)).area,...
            '-' SuperRat(ses).units(objectcells(i)).type num2str(objectcells(i)) nameappend]),'tif');
        close(gcf);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now aggregation of cells:
%{
Object cells: has to fire at least one on at least half of the sampling
events and has to be bootstrap validated different firing rates across
odors

sampling pool: odor selective and has to fire for at least one run?

place cells: has to have a place field in at least one run
    A place field means cell fires at least 3 hz and firing field as at
    most half of the trajectory (top 75% of the rate bins)

1. do odor selective cells show different spatial characteristics than non
selective cells? e.g.
    a. number of place fields
    b. spatial info
    c. mean rate
    d. max traj rate
    e. sparsity
    f. trajectory selectivity
    g.

2. does the odor selectivity relate to the trajectory selectivity?
    e.g. odor selectivity score vs traj selectivty score (diff over sum, grand mean)
	or dprime for the trial rate differences.
    This can be absolute selectivity, or trajectory signed (L odor with L
    traj?



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

%% This runs summary statistics on the interaction between spatial selectivity
% and object selectivity.  Run this separately for each region (PFC and
% HPC)

location='CA1';
%location='PFC';

% whose lin place field plots have information?
% lets just aggregate all the odor selective cells and get the % of them
% who have spatial information, and the % of them who have place fields

allstats=[];

% gather all the pyramidal cells
ispyram=cell2mat(cellfun(@(a) contains(a,'pyr'), {SuperUnits.type}, 'UniformOutput', false)');

% has to have a max rate on any trajectory of at least 2 hz
ActivePop=cell2mat(cellfun(@(a) max(a.PFmax,[],2)>2, {SuperUnits.FieldProps}, 'UniformOutput', false)');

%ActivePop=cell2mat(cellfun(@(a) any(a(:)), {SuperUnits.FiresDuringRun},'UniformOutput', false));
isHere=cell2mat(cellfun(@(a) contains(a,location), {SuperUnits.area},'UniformOutput', false));




usecell=ActivePop(:) & isHere(:) & ispyram(:);

% 1. which are obj cells
allstats=cell2mat(cellfun(@(a) a(:,3), {SuperUnits(usecell).OdorSelective}, 'UniformOutput', false)');
% 2. which have place fields?
allstats(:,2)=cell2mat(cellfun(@(a) sum(a(1:2)), {SuperUnits(usecell).PFexist}, 'UniformOutput', false)');
% howabout spatial information and p values?
% 3. Mean rate during runs
allstats(:,3)=cell2mat(cellfun(@(a) mean(a{1}(:)), {SuperUnits(usecell).LinPlaceFields}, 'UniformOutput', false)');
% 4. information in bits
allstats(:,4)=cell2mat(cellfun(@(a) max(a.info,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 5. sparsity
allstats(:,5)=cell2mat(cellfun(@(a) max(a.sparsity,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 6. place field size
allstats(:,6)=cell2mat(cellfun(@(a) min(a.PFsize,[],2), {SuperUnits(usecell).FieldProps}, 'UniformOutput', false)');
% 7. how about the peak rate
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
% 14 is whether the peak is high enough (i'd say both have to have peaks
% above like 2 hz?
allstats(:,14)=cellfun(@(a) max(a{1}(1,:))>3 && max(a{1}(2,:))>3, {SuperUnits(usecell).LinPlaceFields});
% 15 is the correlation between trajectories 1 and 2
allstats(:,15)=cell2mat(cellfun(@(a) corr(a{1}(1,:)',a{1}(2,:)'), {SuperUnits(usecell).LinPlac
    eFields}, 'UniformOutput', false))';
% 16 is does it fire on each trajectory? (already has to fire +2 on a% ma
% trajectory
allstats(:,16)=cell2mat(cellfun(@(a) sum(a(1:2)), {SuperUnits(usecell).FiresDuringRun}, 'UniformOutput', false)');


%%

% cut out all nanned out sessions
nancells=isnan(allstats(:,1));
allstats(nancells,:)=[];

% howabout cap mean rate overall

% cant be ins and have to spike af ew times?
pcells=allstats(:,3)<8; % cant be INs
%pcells=allstats(:,2)==1; %
%pcells=ones(size(allstats,1),1)==1;


% now plot
figure('Position',[150 150 1000 1000]);
xnamecell={'Not','Odor selective'};
xnames=categorical(xnamecell); xnames=reordercats(xnames,xnamecell);

% %% of cells w fields
subplot(3,3,1);
hb=bar(xnames,100*[nanmean(allstats(allstats(:,1)==0,2)==1) nanmean(allstats(allstats(:,1)==1,2)==1);...
    nanmean(allstats(allstats(:,1)==0,2)==2) nanmean(allstats(allstats(:,1)==1,2)==2)],'stacked');
ylabel('% of cells with a place field'); title(sprintf('Non Odor cells %d / %d, Odor %d/%d',...
    nansum(allstats(allstats(:,1)==0,2)>0), nansum(allstats(:,1)==0),...
    nansum(allstats(allstats(:,1)==1,2)>0), nansum(allstats(:,1)==1)));
legend('One pf','Both PF'); set(gca,'ylim',[0 100]);

subplot(3,3,2);
% mean rate
boxScatterplot(allstats(pcells,3),allstats(pcells,1),'position',[],'xLabels',xnamecell,...
    'Ylabel','Overall Mean rate (Hz)');
title(sprintf('Mean Rate: %.2f & %.2f \n ranksum p = %.4f',...
    nanmean(allstats(pcells & allstats(:,1)==0,3)),nanmean(allstats(pcells & allstats(:,1)==1,3)),...
    ranksum(allstats(pcells & allstats(:,1)==1,3),allstats(pcells & allstats(:,1)==0,3))));


subplot(3,3,3);
% seventh is peak rate
boxScatterplot(allstats(pcells,7),allstats(pcells,1),'position',[],'xLabels',xnamecell,...
    'Ylabel','PF peak rate');
title(sprintf('Peak Rate: %.2f & %.2f \n ranksum p = %.4f',...
    nanmean(allstats(pcells & allstats(:,1)==0,7)),nanmean(allstats(pcells & allstats(:,1)==1,7)),...
    ranksum(allstats(pcells & allstats(:,1)==1,7),allstats(pcells & allstats(:,1)==0,7))));

%%%%% only use cells w place fields from here on out
% mean spatial information
subplot(3,3,4);
boxScatterplot(allstats(pcells,4),allstats(pcells,1),'position',[],'xLabels',xnamecell,...
    'Ylabel','bits per spike');
title(sprintf('Information in bits: %.2f & %.2f \n ranksum p=%.4f',...
    nanmean(allstats(pcells & allstats(:,1)==0,4)),nanmean(allstats(pcells & allstats(:,1)==1,4)),...
    ranksum(allstats((allstats(:,1)==1 & pcells),4),allstats((allstats(:,1)==0 & pcells),4))));




%{
subplot(3,3,5);
% fifth col is max sparsity
boxScatterplot(log2(allstats(pcells,5)),allstats(pcells,1),'position',[],'xLabels',xnamecell,...
    'Ylabel','log sparsity');
title(sprintf('Sparsity: %.2f & %.2f \n ranksum p = %.4f',...
    nanmean(allstats(pcells & allstats(:,1)==0,5)),nanmean(allstats(pcells & allstats(:,1)==1,5)),...
    ranksum(allstats(pcells & allstats(:,1)==1,5),allstats(pcells & allstats(:,1)==0,5))));

box off;
    %}
    
    subplot(3,3,5);
    % howabout max/mean
    maxmean=allstats(:,7)-allstats(:,3);
    boxScatterplot(maxmean(pcells),allstats(pcells,1),'position',[],'xLabels',xnamecell,...
        'Ylabel','log sparsity');
    title(sprintf('Max over mean %.2f & %.2f \n ranksum p = %.4f',...
        nanmean(maxmean(pcells & allstats(:,1)==0)),nanmean(maxmean(pcells & allstats(:,1)==1)),...
        ranksum(maxmean(pcells & allstats(:,1)==1),maxmean(pcells & allstats(:,1)==0))));
    
    
    
    subplot(3,3,6);
    % Sixth is pfsize
    % but remove the cells whose peak rate is below 1
    notnans=allstats(:,7)>=2 & allstats(:,6)<99 & pcells;
    %okcells=pcells;
    boxScatterplot(allstats(notnans,6),allstats(notnans,1),'position',[],'xLabels',xnamecell,...
        'Ylabel','PF sizes');
    title(sprintf('Field Size: %.2f & %.2f \n ranksum p = %.4f',...
        nanmean(allstats(pcells & allstats(:,1)==0,6)),nanmean(allstats(pcells & allstats(:,1)==1,6)),...
        ranksum(allstats(allstats(notnans,1)==1,6),allstats(allstats(notnans,1)==0,6))));
    
    
    
    subplot(3,3,7);
    % eighth trajectory score (Abs val)
    boxScatterplot(abs(allstats(:,8)),allstats(:,1),'position',[],'xLabels',xnamecell,...
        'Ylabel','Run Selectivity Score');
    title(sprintf('Path Selectivity: %.4f & %.4f \n ranksum p = %.6f',...
        nanmean(abs(allstats(pcells & allstats(:,1)==0,8))),nanmean(abs(allstats(pcells & allstats(:,1)==1,8))),...
        ranksum(abs(allstats(allstats(:,1)==1,8)),abs(allstats(allstats(:,1)==0,8)))));
    
    if contains(location,'CA1')
        sgtitle({'Hippocampal Pyramidal cells'; 'Object Selective Cells vs non selective'; ''});
    else
        sgtitle({'Prefrontal Pyramidal cells'; 'Object Selective Cells vs non selective'; ''});
    end
    
    
    %{
% for ca1, it looks like object selective cells are more likely to be path
% selective, lets see the 2x2 square
% so lets see if theres an interaction between which path and which odor
%    1 is pval odor       10 is dir odor   11 is pval traj 8 is traj dir
odorpr=nanmean(allstats(:,1)==1 & allstats(:,10)>0 & allstats(:,11)>0 & allstats(:,8)>0); % odor 1 selective %
odorchi=sum(allstats(:,1)==1 & allstats(:,10)>0 & allstats(:,11)>0 & allstats(:,8)>0);
odorchi(1,2)=sum(allstats(:,1)==1 & allstats(:,10)>0 & allstats(:,11)>0 & allstats(:,8)<0);
odorchi(2,1)=sum(allstats(:,1)==1 & allstats(:,10)<0 & allstats(:,11)>0 & allstats(:,8)>0);
odorchi(2,2)=sum(allstats(:,1)==1 & allstats(:,10)<0 & allstats(:,11)>0 & allstats(:,8)<0);

% basically what are the overlaps between object selective and trajectory selective?
nanmean(allstats(allstats(:,1)==1,11)<.05)
nanmean(allstats(allstats(:,1)==0,11)<.05)
    %}
    
    % do odor selective cells show trajectory dependent firing elsewhere
    % e.g. is that odor cell firing later with the same selectivity?
    % so odor selectivity will be the dprime or the diff over sum score
    %  ** this needs to be signed by which is larger
    % trajectory splitting will be
    
    
    subplot(3,3,8);
    % Do odor selective units select the same subesequent route??
    
    
    plot(allstats(pcells & allstats(:,1)==0,9),allstats(pcells & allstats(:,1)==0,8),'.');
    hold on;
    plot(allstats(pcells & allstats(:,1)==1,9),allstats(pcells & allstats(:,1)==1,8),'.');
    
    %plot((corstats(corstats(:,2)==0,1)),(corstats(corstats(:,2)==0,3)),'.');
    %hold on;
    %plot((corstats(corstats(:,2)==1,1)),(corstats(corstats(:,2)==1,3)),'.');
    %legend('Not','Odor Selective')
    
    plot(get(gca,'XLim'),[0 0],'r');
    hold on; plot([0 0],get(gca,'YLim'),'r'); axis tight;
    ylabel('Right traj pref <  > Left traj pref');
    xlabel('Right odor pref <  > Left odor Pref');
    notnans=~isnan(allstats(pcells & allstats(:,1)==0,9)) & ~isnan(allstats(pcells & allstats(:,1)==0,8));
    %[a,b]=corr((allstats(okcells & allstats(:,2)==1,1)),(allstats(okcells & allstats(:,2)==1,3)));
    [a,b]=corr(allstats(pcells & allstats(:,1)==0,9),allstats(pcells & allstats(:,1)==0,8));
    title(sprintf('Directional selectivity rho=%.2f \n Spearman P= %.4f',a,b));
    % looks like an inverse relationship...
    % which would make sense if the odor selective cells have worse place
    % fields
    usecorr=1;
    if ~usecorr
        subplot(3,3,9);
        plot(allstats(pcells & allstats(:,1)==0,9),allstats(pcells & allstats(:,1)==0,8),'.');
        hold on;
        plot(allstats(pcells & allstats(:,1)==1,9),allstats(pcells & allstats(:,1)==1,8),'.');
        
        axis tight; ylabel('Trajectory score'); xlabel('Odor Preference (Log_1_0)');
        notnans=~isnan(allstats(pcells & allstats(:,1)==0,9)) &...
            ~isnan(allstats(pcells & allstats(:,1)==0,8));
        [a,b]=corr(abs(allstats(pcells & allstats(:,1)==0,9)),abs(allstats(pcells & allstats(:,1)==0,8)));
        title(sprintf('Nondirectional selectivity rho=%.2f \n Spearman P= %.4f',a,b));
    else
        subplot(3,3,9);
        plot(log10(abs(allstats(pcells & allstats(:,1)==0,9))),allstats(pcells & allstats(:,1)==0,15),'.');
        
        hold on;
        plot(log10(abs(allstats(pcells & allstats(:,1)==1,9))),allstats(pcells & allstats(:,1)==1,15),'.');
        
        plot(get(gca,'XLim'),[0 0],'r');
        hold on; plot([0 0],get(gca,'YLim'),'r'); axis tight;
        
        ylabel('Trajectory Correlation');
        xlabel('Odor Preference (Log_1_0)');
        notnans=~isnan(allstats(:,9)) & ~isnan(allstats(:,15));
        [a,b]=corr(abs(allstats(pcells & notnans,9)),allstats(pcells & notnans,15));
        title(sprintf('Path similarity vs odor selectivity rho=%.2f \n Spearman P= %.4f',a,b));
        
    end
    
    %% is this effect in sp 9 caused by firing rate?
    % it doenst look to be so but it definitely could be
    figure;
    scatter(log10(abs(allstats(pcells,9))),allstats(pcells,15),25,rescale(log(allstats(pcells,3))),'filled');
    
    hold on; plot([0 0],get(gca,'YLim'),'r'); axis tight;
    
    ylabel('Trajectory Correlation');
    xlabel('Odor Preference (Log_1_0)');
    notnans=~isnan(allstats(:,9)) & ~isnan(allstats(:,15));
    [a,b]=corr(abs(allstats(pcells & notnans,9)),allstats(pcells & notnans,15));
    title(sprintf('Path similarity vs odor selectivity rho=%.2f \n Spearman P= %.4f',a,b));
    
    % doesnt look like it
    
    %% where do the place fields of odor selective cells tend to exist?
    % take all cells with place fields, normalize by max, concatenate and show
    % pf density.
    % then split those cells into odor selelctive cells and non odor selective
    
    
    % so it looks to me that the odor selective cells have a higher run
    % selectivity, but the spatial organization of the spikes are more
    % correlated.
    
    % So the idea might be that odor selective cells have a similar spatial
    % layout, but their rates are more reliably different.  Almost like its a
    % rate remapping effect over a global remapping effect
    
    % so one method is from the moser paper: Liu et al 2013
    % they used three methods to compare conditions:
    
    % 1. change in rate
    % 2. spatial correlation
    % 3. PV correlation
    % the gist is that the pv correlation changed a lot
    
    % so the possible measures that will flush this out
    % peak shift btwn trajs
    % peak difference btwn trajs
    
    % need to figure out how they calc pv correlation
    % pv correlation was done by catting the mean place fields across x positions
    % and then correlating the y stacks so you'll have a correlation across the
    % x that you can plot
    
    figure;
    subplot(3,3,[4 7]);
    % 12 is difference in peak rates
    % ok cells are those who have peaks to compare and are pyrams
    notnans=allstats(:,14)==1 & pcells;
    boxScatterplot(abs(allstats(notnans,12)),allstats(notnans,1),'position',[],'xLabels',xnamecell,...
        'Ylabel','Difference between pf peak rates');
    title(sprintf('Diff in peak rates: %.4f & %.4f \n ranksum p = %.6f',...
        nanmean(abs(allstats(allstats(notnans,1)==0,12))),nanmean(abs(allstats(allstats(notnans,1)==1,12))),...
        ranksum(abs(allstats(allstats(:,1)==1,12)),abs(allstats(allstats(:,1)==0,12)))));
    
    subplot(3,3,[5 8]);
    notnans=allstats(:,14)==1 & pcells;
    boxScatterplot(abs(allstats(notnans,13)),allstats(notnans,1),'position',[],'xLabels',xnamecell,...
        'Ylabel','Distance btwn Peaks');
    title(sprintf('Distance btwn peaks: %.4f & %.4f \n ranksum p = %.6f',...
        nanmean(abs(allstats(allstats(notnans,1)==0,13))),nanmean(abs(allstats(allstats(notnans,1)==1,13))),...
        ranksum(abs(allstats(allstats(:,1)==1,13)),abs(allstats(allstats(:,1)==0,13)))));
    
    subplot(3,3,[6 9]);
    notnans=allstats(:,14)==1 & pcells;
    boxScatterplot(abs(allstats(notnans,8)),allstats(notnans,1),'position',[],'xLabels',xnamecell,...
        'Ylabel','Difference in mean traj rates');
    title(sprintf('Difference in mean traj rates: %.4f & %.4f \n ranksum p = %.6f',...
        nanmean(abs(allstats(allstats(notnans,1)==0,8))),nanmean(abs(allstats(allstats(notnans,1)==1,8))),...
        ranksum(abs(allstats(allstats(:,1)==1,8)),abs(allstats(allstats(:,1)==0,8)))));
    
    
    
    
    if contains(location,'CA1')
        sgtitle({'Hippocampal Pyramidal cells: Object Selective Cells vs non selective'});
    else
        sgtitle({'Prefrontal Pyramidal cells: Object Selective Cells vs non selective'});
    end
    
    
    %% Liu Pop vector encoding
    
    % concatenate all the trajectories for left, and separately for the right
    % z-score the rate vectors, then correlate left vs right
    % this will yield a correlation vs position for ensembles that are
    % selective, and that for ensembles that arent.
    
    pvleft=zscore(cell2mat(cellfun(@(a) a{1}(1,:), {SuperUnits(usecell).LinPlaceFields}, 'UniformOutput', false)'),1,2);
    pvright=zscore(cell2mat(cellfun(@(a) a{1}(2,:), {SuperUnits(usecell).LinPlaceFields}, 'UniformOutput', false)'),1,2);
    pvleft(nancells,:)=[];
    pvright(nancells,:)=[];
   
   % this gets the pvleft and pvright to the same place as the allstats
   
   
   
    % gotta remove nans, and parse this region, can also only use pyrams if you
    % want
    notnans=sum(isnan(pvleft),2)==0 & sum(isnan(pvright),2)==0;
    isobj=
    
    usecell=notnans & myregion';
    
    cormat=corr(pvleft(usecell & isobj,:),pvright(usecell & isobj,:));
    cormat2=corr(pvleft(usecell & ~isobj,:),pvright(usecell & ~isobj,:));
    
    
    figure; subplot(1,2,1);
    imagesc(cormat); title('object selective cells');
    subplot(1,2,2); imagesc(cormat2); title('non object selective cells');
    if contains(location,'CA1')
        sgtitle({'Hippocampal Pyramidal cells: Object Selective Cells vs non selective'});
    else
        sgtitle({'Prefrontal Pyramidal cells: Object Selective Cells vs non selective'});
    end
    
    ranksum(cormat(eye(100)==1),cormat2(eye(100)==1));
    
    % I think the real control is to downsample the first cormat and see if it
    % works...
    
