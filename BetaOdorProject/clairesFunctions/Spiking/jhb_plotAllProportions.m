%%jhb_plotAllProportions

% this function will first plot the phase locking proportions for local
% coherence, then in a stacked bar for cross-regional coherence


% so the variables I need:
%1. npCells_CA1/PFC old and npInt_CA1/PFC old
%2. under analysisacrossanimals, phaselocking, plcells_beta and rr, all
%regions
%3. pyrs vs ins

[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

%% grab our underlying data
if ~exist('cellpops_CA1','var')
    clear
    topDir = cs_setPaths();
    regions = {'CA1','PFC'};

    
    for r=1:length(regions)
        load([topDir,'AnalysesAcrossAnimals\cellPopulations_',regions{r}],'cellpops');
        eval(['cellpops_',regions{r}, '=cellpops; openvar(''cellpops_', regions{r}, '''); clear cellpops']);
    end
end

%% Proportionof Task Responsive cells Locked to Any regions Rhythm

% a paired bar is n by 2
bands={'beta','resp'};
celltypes={'pyr','int'};

cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];
bandcolors=[rgbcolormap('MidNightBlue'); rgbcolormap('Magenta')];


barProps=nan(4,2);
barCts=nan(4,2);
npCts=nan(4,2);
for barc=1:length(bands)
    % first get ca1 pyr then int
    for region=1:length(regions) % 12, then 34
        
        % first pyrams
        eval(['npCells=cellpops_' regions{region} '.npPyrams;']);
        eval(['cohCells=npCells(ismember(cellpops_' regions{region} '.' bands{barc} 'PhaseLockedAny,npCells,''rows''),:);']);
        
        npCts((region-1)*2+1,barc)=size(npCells,1);
        barCts((region-1)*2+1,barc)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' regions{region} '.npInt;']);
        eval(['cohInt=npInt(ismember(cellpops_' regions{region} '.' bands{barc} 'IntPhaseLockedAny,npInt,''rows''),:);']);
        
        npCts(region*2,barc)=size(npInt,1);
        barCts(region*2,barc)=size(cohInt,1);
        
    end
end

barProps=barCts./npCts;
figure('Position',[1578, 853, 313, 316]);
hb=bar(barProps,1); box off;
ylabel(sprintf('Proportion of \n Task Responsive Cells'));
set(gca,'XTick',1:4,'XTickLabel',{'Pyr','Int'});
hold on; plot([0.5 4.5],[.05 .05],'k--'); xlim([.5 4.5]);
xlabel('CA1              PFC'); legend(hb,bands); legend('boxoff')
for barc=1:2
    set(hb(barc),'LineStyle','none','FaceColor',bandcolors(barc,:));
end
title('Locked to Any Region')
proptable=array2table([npCts(:,1) barCts],'VariableNames',{'taskResponsive','BetaCoh','RespCoh'},...
    'RowNames',{'CA1Pyr','CA1Int','PFCPyr','PFCInt'});
openvar('proptable');


%% Proportion of Task Responsive cells Locked to LOCAL rhythm
% right now unused
% a paired bar is n by 2
bands={'beta','resp'};
celltypes={'pyr','int'};

cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];
bandcolors=[rgbcolormap('MidNightBlue'); rgbcolormap('Magenta')];


barProps=nan(4,2);
barCts=nan(4,2);
npCts=nan(4,2);
for barc=1:length(bands)
    % first get ca1 pyr then int
    for region=1:length(regions) % 12, then 34
        
        % first pyrams
        eval(['npCells=cellpops_' regions{region} '.npPyrams;']);
        eval(['cohCells=npCells(ismember(cellpops_' regions{region} '.' bands{barc} 'PhaseLockedTo' regions{region} ',npCells,''rows''),:);']);
        
        npCts((region-1)*2+1,barc)=size(npCells,1);
        barCts((region-1)*2+1,barc)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' regions{region} '.npInt;']);
        eval(['cohInt=npInt(ismember(cellpops_' regions{region} '.' bands{barc} 'IntPhaseLockedTo' regions{region} ',npInt,''rows''),:);']);
        
        npCts(region*2,barc)=size(npInt,1);
        barCts(region*2,barc)=size(cohInt,1);
        
    end
end

barProps=barCts./npCts;
figure('Position',[1578, 853, 313, 316]);
hb=bar(barProps,1); box off;
ylabel(sprintf('Proportion of \n Task Responsive Cells'));
set(gca,'XTick',1:4,'XTickLabel',{'Pyr','Int'},'YTick',0:.2:1);
hold on; plot([0.5 4.5],[.05 .05],'k--'); xlim([.5 4.5]);
xlabel('CA1              PFC'); legend(hb,bands); legend('boxoff')
for barc=1:2
    set(hb(barc),'LineStyle','none','FaceColor',bandcolors(barc,:));
end
title('Locked to Local LFP')

proptable=array2table([npCts(:,1) barCts],'VariableNames',{'taskResponsive','BetaCoh','RespCoh'},...
    'RowNames',{'CA1Pyr','CA1Int','PFCPyr','PFCInt'});
openvar('proptable');

%% Beta coherence pyr in both cell regions to all LFPregions
% need to find a way to get these nonoverlapping... cells are counted twice
% here

% stacked bar rows are x pos, columns are color series in bar

cellcolors=[rgbcolormap('Salmon'); rgbcolormap('DarkTurquoise'); rgbcolormap('DarkOrange')];
bandcolors=[rgbcolormap('MidNightBlue'); rgbcolormap('Magenta')];

cellregions={'CA1','PFC'};
eegregions={'CA1','PFC','OB'};
freq='resp';


barProps=nan(4,3); % rows are ca1pyr ca1 in pfc pyr.... cols are ca1 pfc ob
barCts=nan(4,3);
npCts=nan(4,3);
plcells=repmat({[],[],[]},4,1);
for cr=1:length(cellregions)
    
    for er=1:length(eegregions) % 12, then 34
        
        % first pyrams
        eval(['npCells=cellpops_' cellregions{cr} '.npPyrams;']);
        % get all cell info with the pl p value
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            
            load([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',cellregions{cr},'-',eegregions{er},'.mat'])
            cellfilter = '($prayl < 0.05) && ($Nspikes > 0)';
            cells = evaluatefilter(phaselock,cellfilter);
            
            if ~isempty(cells)
                % session, no epoch, tet, unit
                noeps = cells(:,[1 3 4]);
                cells = unique(noeps,'rows');
                % now add animal
                cells = [repmat(a, size(cells,1),1),cells];
                for i=1:size(cells,1)
                    cells(i,5)=er; cells(i,6)=phaselock{cells(i,2)}{1}{cells(i,3)}{cells(i,4)}.prayl;
                end
                plcells{(cr-1)*2+1,er} = [plcells{(cr-1)*2+1,er}; cells];
            end
        end
        
        npCts((cr-1)*2+1,er)=size(npCells,1);
        %barCts((cr-1)*2+1,er)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' cellregions{cr} '.npInt;']);
       for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            
            load([animDir,'PhaseLocking\Interneurons\',animal,'phaselock_',freq,'_',cellregions{cr},'-',eegregions{er},'.mat'])
            cellfilter = '($prayl < 0.05) && ($Nspikes > 0)';
            cells = evaluatefilter(phaselock,cellfilter);
            
            if ~isempty(cells)
                % session, no epoch, tet, unit
                noeps = cells(:,[1 3 4]);
                cells = unique(noeps,'rows');
                % now add animal
                cells = [repmat(a, size(cells,1),1),cells];
                for i=1:size(cells,1)
                    cells(i,5)=er; cells(i,6)=phaselock{cells(i,2)}{1}{cells(i,3)}{cells(i,4)}.prayl;
                end
                plcells{cr*2,er} = [plcells{cr*2,er}; cells];
            end
        end
        npCts(cr*2,er)=size(npInt,1);
        %barCts(cr*2,er)=size(cohInt,1);
        
    end
end
% now generate bar counts from plcells

% now reform our plcells rows so that each cell only occurs once
% the cells are 1. rat 2. session 3. tet, 4. unit, 5. target region 6.
% prayl
% so for each row
for ct=1:4 % celltype
    allplcells=cell2mat(plcells(ct,:)');
    [unitUnique,~,occurrences]=unique(allplcells(:,1:4),'rows');
    unitUnique(:,5:6)=nan;
    for k=1:size(unitUnique,1)
        cans=allplcells(occurrences==k,:);
        [~,winner]=min(cans(:,6));
        unitUnique(k,:)=cans(winner,:);
    end
    
    barCts(ct,:)=accumarray(unitUnique(:,5),1,[3,1]);
end

barProps=barCts./npCts;
figure('Position',[1578, 853, 313, 316]);
hb=bar(barProps,'stacked'); box off;
ylabel(sprintf('Proportion of \n Task Responsive Cells'));
set(gca,'XTick',1:4,'XTickLabel',{'Pyr','Int'},'YTick',0:.2:1);
hold on; plot([0.5 4.5],[.05 .05],'k--'); xlim([.5 4.5]);

xlabel('CA1              PFC'); 
legend(hb,eegregions); legend('boxoff')
for barc=1:3
    set(hb(barc),'FaceColor',cellcolors(barc,:));
end
title(sprintf('%s Locking',freq))

proptable=array2table([npCts(:,1) barCts],'VariableNames',{'taskResponsive',['CA1' freq],['PFC' freq],['OB' freq]},...
    'RowNames',{'CA1Pyr','CA1Int','PFCPyr','PFCInt'});
openvar('proptable');

%% Resp coherence pyr in both cell regions to all LFPregions

% stacked bar rows are x pos, columns are color series in bar

cellcolors=[rgbcolormap('Salmon'); rgbcolormap('DarkTurquoise'); rgbcolormap('DarkOrange')];
bandcolors=[rgbcolormap('MidNightBlue'); rgbcolormap('Magenta')];

cellregions={'CA1','PFC'};
eegregions={'CA1','PFC','OB'};

barProps=nan(4,3); % rows are ca1pyr ca1 in pfc pyr.... cols are ca1 pfc ob
barCts=nan(4,3);
npCts=nan(4,3);
for cr=1:length(cellregions)
    
    for er=1:length(eegregions) % 12, then 34
        
        % first pyrams
        eval(['npCells=cellpops_' cellregions{cr} '.npPyrams;']);
        eval(['cohCells=npCells(ismember(cellpops_' cellregions{cr} '.respPhaseLockedTo' eegregions{er} ',npCells,''rows''),:);']);
        
        npCts((cr-1)*2+1,er)=size(npCells,1);
        barCts((cr-1)*2+1,er)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' cellregions{cr} '.npInt;']);
        try
            eval(['cohInt=npInt(ismember(cellpops_' cellregions{cr} '.respIntPhaseLockedTo' eegregions{er} ',npInt,''rows''),:);']);
        catch
            cohInt=[]; % this happens if there are no npInt cells, and the ismember fx cant work
        end
        npCts(cr*2,er)=size(npInt,1);
        barCts(cr*2,er)=size(cohInt,1);
        
    end
end

barProps=barCts./npCts;
figure('Position',[1578, 853, 313, 316]);
hb=bar(barProps,'stacked'); box off;
ylabel(sprintf('Proportion of \n Task Responsive Cells'));
set(gca,'XTick',1:4,'XTickLabel',{'Pyr','Int'},'YTick',0:.2:1);
hold on; plot([0.5 4.5],[.05 .05],'k--'); xlim([.5 4.5]);

xlabel('CA1              PFC'); 
legend(hb,eegregions); legend('boxoff')
for barc=1:3
    set(hb(barc),'FaceColor',cellcolors(barc,:));
end
title('Beta Locking')
%% now overlap between beta and rr locking (any) for pyramidal cells
       
bandcolors=[rgbcolormap('MidNightBlue'); rgbcolormap('Magenta')];
bands={'beta','resp'};
% pyrams first
figure;

barProps=nan(8,1); % rows are ca1pyr ca1 in pfc pyr.... cols are ca1 pfc ob
barCts=nan(8,1);
npCts=nan(2,1);
% so the bar is ca1: beta rr both chance
for cr=1:length(cellregions)
    eval(['npCells=cellpops_' cellregions{cr} '.npPyrams;']);
    npCts(cr)=size(npCells,1);
    cohBoth={};
    for band=1:length(bands)
        eval(['cohCells=npCells(ismember(cellpops_' cellregions{cr} '.' bands{band},...
            'PhaseLockedAny,npCells,''rows''),:);']);
        barCts(band+(cr-1)*4)=size(cohCells,1);
        cohBoth{band}=cohCells;
    end
    dualCoh=intersect(cohBoth{1},cohBoth{2},'rows');
    barCts((cr-1)*4+3)=size(dualCoh,1);
    nullProb=barCts(1:2+(cr-1)*4)
    [~,CI]=binofit(,npCts(cr));
    




