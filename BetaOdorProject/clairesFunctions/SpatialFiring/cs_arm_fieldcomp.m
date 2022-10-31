% place fields on center stem only- comparison between this selectivity vs odor selectivity
% loop through all place cells (load linfields). find ones that have field on center stem on outbound trials.
% use linpos to determine where the center stem ends? 
% for cells that have fields on center stem, save in new struct, splitter cells
% 
% get the field SI, and the odor SI. 
% calculate CC and significance 


[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
dataDir = [topDir,'AnalysesAcrossAnimals\'];
region = 'CA1';
load([dataDir,'placeCells_',region]);
load([dataDir, 'placemaps']);
comps = [];
for a = [1 2 3 4 8]
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    animcells = placecells(placecells(:,1)==a,(2:4));
    
    days = unique(animcells(:,1));
    for d = days'
        linfields = loaddatastruct(animDir, animal, 'linfields',d);
        linpos = loaddatastruct(animDir, animal, 'linpos',d);
        cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
        
        epochs = cs_getRunEpochs(animDir, animal, 'odorplace',d);
        epochs = epochs(:,2);
        cells = animcells(animcells(:,1)==d,(2:3));
        
        for c = 1:size(cells,1)
            cell = cells(c,:);
            cellepochs = cs_findGoodEpochs(cellinfo{d}, 'SI',cell);
            if isempty(cellepochs) %cell did not spike during nosepoke
                continue
            end
            
            if ~isfield(cellinfo{d}{cellepochs(1)}{cell(1)}{cell(2)},'SI')
                continue
            else
            odorSI = cellinfo{d}{cellepochs(1)}{cell(1)}{cell(2)}.SI;
            end
            
            fields_l = [];
            fields_r = [];
            
            for e = epochs'
                if length(linfields{d}{e}{cell(1)})>=cell(2) && ~isempty(linfields{d}{e}{cell(1)}{cell(2)})
                stemlength = linpos{d}{e}.segmentInfo.segmentLength(1);
                 armfield_l = linfields{d}{e}{cell(1)}{cell(2)}{1}(linfields{d}{e}{cell(1)}{cell(2)}{1}(:,1) >= stemlength,5);
                 armfield_r = linfields{d}{e}{cell(1)}{cell(2)}{3}(linfields{d}{e}{cell(1)}{cell(2)}{3}(:,1) >= stemlength,5);
                fields_l = stack(fields_l,armfield_l');
                fields_r = stack(fields_r, armfield_r');
                else
                    continue
                end
            end
            
            field_l = nanmean(fields_l,1);
            field_r = nanmean(fields_r,1);
            
            fields = [field_l,field_r];
            [peak,peakind] = max(fields);
            if peak >= 3 %cell has field on stem
                meanfr_l = nanmean(field_l);
                meanfr_r = nanmean(field_r);
                trajSI = (meanfr_l - meanfr_r) / (meanfr_l + meanfr_r);
                comps = [comps; odorSI, trajSI];
                
                
%                 cellnum = find(ismember(placecells, [a, d, cell],'rows'));
%                 field = placemaps{cellnum};
%                 imagesc(field)
%                 colorbar
%                 title(num2str(odorSI));
%                 pause
                
            end
            
            
            
        end
    end
end

%Also plot individual cells- may need to limit to only stem, not including
%turn. 

[bad,~] = find(isnan(comps));
    comps(bad,:) = [];
    
    figure,
    plot(comps(:,1), comps(:,2), 'k.');
    hold on
    fit = polyfit(comps(:,1), comps(:,2),1);
    plot([-1 1], polyval(fit,[min(comps(:,1)), max(comps(:,1))]))
    axis([-1 1 -1 1]);
    
    title(region)
    
    [CC,p] = corrcoef(comps(:,1), comps(:,2));
    R = CC(1,2)
    p = p(1,2)
    text(0.5,0.7,{['R = ' num2str(R)],;['p = ' num2str(p)]})
    figfile = [figDir, 'PlaceFields\SI-ArmCorrelation_', region];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);