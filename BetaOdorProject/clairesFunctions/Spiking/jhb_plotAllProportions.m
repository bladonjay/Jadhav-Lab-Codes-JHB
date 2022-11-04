%%jhb_plotAllProportions

% this function will first plot the phase locking proportions for local
% coherence, then in a stacked bar for cross-regional coherence


% so the variables I need:
%1. npCells_CA1/PFC old and npInt_CA1/PFC old
%2. under analysisacrossanimals, phaselocking, plcells_beta and rr, all
%regions
%3. pyrs vs ins

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
for band=1:length(bands)
    % first get ca1 pyr then int
    for region=1:length(regions) % 12, then 34
        
        % first pyrams
        eval(['npCells=cellpops_' regions{region} '.npPyrams;']);
        eval(['cohCells=npCells(ismember(cellpops_' regions{region} '.' bands{band} 'PhaseLockedAny,npCells,''rows''),:);']);
        
        npCts((region-1)*2+1,band)=size(npCells,1);
        barCts((region-1)*2+1,band)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' regions{region} '.npInt;']);
        eval(['cohInt=npInt(ismember(cellpops_' regions{region} '.' bands{band} 'IntPhaseLockedAny,npInt,''rows''),:);']);
        
        npCts(region*2,band)=size(npInt,1);
        barCts(region*2,band)=size(cohInt,1);
        
    end
end

barProps=barCts./npCts;
figure('Position',[1578, 853, 313, 316]);
hb=bar(barProps,1); box off;
ylabel(sprintf('Proportion of \n Task Responsive Cells'));
set(gca,'XTick',1:4,'XTickLabel',{'Pyr','Int'});
hold on; plot([0.5 4.5],[.05 .05],'k--'); xlim([.5 4.5]);
xlabel('CA1              PFC'); legend(hb,bands); legend('boxoff')
for band=1:2
    set(hb(band),'LineStyle','none','FaceColor',bandcolors(band,:));
end
title('Locked to Any Region')


%% Proportion of Task Responsive cells Locked to LOCAL rhythm

% a paired bar is n by 2
bands={'beta','resp'};
celltypes={'pyr','int'};

cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];
bandcolors=[rgbcolormap('MidNightBlue'); rgbcolormap('Magenta')];


barProps=nan(4,2);
barCts=nan(4,2);
npCts=nan(4,2);
for band=1:length(bands)
    % first get ca1 pyr then int
    for region=1:length(regions) % 12, then 34
        
        % first pyrams
        eval(['npCells=cellpops_' regions{region} '.npPyrams;']);
        eval(['cohCells=npCells(ismember(cellpops_' regions{region} '.' bands{band} 'PhaseLockedTo' regions{region} ',npCells,''rows''),:);']);
        
        npCts((region-1)*2+1,band)=size(npCells,1);
        barCts((region-1)*2+1,band)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' regions{region} '.npInt;']);
        eval(['cohInt=npInt(ismember(cellpops_' regions{region} '.' bands{band} 'IntPhaseLockedTo' regions{region} ',npInt,''rows''),:);']);
        
        npCts(region*2,band)=size(npInt,1);
        barCts(region*2,band)=size(cohInt,1);
        
    end
end

barProps=barCts./npCts;
figure('Position',[1578, 853, 313, 316]);
hb=bar(barProps,1); box off;
ylabel(sprintf('Proportion of \n Task Responsive Cells'));
set(gca,'XTick',1:4,'XTickLabel',{'Pyr','Int'},'YTick',0:.2:1);
hold on; plot([0.5 4.5],[.05 .05],'k--'); xlim([.5 4.5]);
xlabel('CA1              PFC'); legend(hb,bands); legend('boxoff')
for band=1:2
    set(hb(band),'LineStyle','none','FaceColor',bandcolors(band,:));
end
title('Locked to Local LFP')
%% Beta coherence pyr in both cell regions to all LFPregions

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
        eval(['cohCells=npCells(ismember(cellpops_' cellregions{cr} '.betaPhaseLockedTo' eegregions{er} ',npCells,''rows''),:);']);
        
        npCts((cr-1)*2+1,er)=size(npCells,1);
        barCts((cr-1)*2+1,er)=size(cohCells,1);
        
        % now ints
        eval(['npInt=cellpops_' cellregions{cr} '.npInt;']);
        try
            eval(['cohInt=npInt(ismember(cellpops_' cellregions{cr} '.betaIntPhaseLockedTo' eegregions{er} ',npInt,''rows''),:);']);
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
for band=1:3
    set(hb(band),'FaceColor',cellcolors(band,:));
end
title('Beta Locking')
