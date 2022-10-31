%cs_placeFields_fieldcomp

%Find field of each place cell (or selective cell) - define region around peak where the FR is greater than threshold, consider this
%the field. Compare FR in this area between R and L trials, get a
%selectivity index for this. Then, compare this to the odor SI for the
%cell, calculate correlation using all cells


clear

thresh = 3; %min FR


animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
[topDir,figDir] = cs_setPaths;

for r = 1:length(regions)
    region = regions{r};
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region,'.mat']);
    %load([topDir,'AnalysesAcrossAnimals\npCells_',region,'.mat']);
    
    fieldcomps = [];
    for a = [1 2 3 4 8] %1:length(animals)  Don't use animals that had truncated track
        animal = animals{a};
        animcells = selectivecells(selectivecells(:,1) == a, 2:4);
        %animcells = npCells(npCells(:,1) == a, 2:4);
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        load([animDir,animal,'cellinfo.mat']);
        days = unique(animcells(:,1));
        
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            load([animDir,animal,'linfields',daystr,'.mat']);
            
            runmatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            
            runmatrix = runmatrix(runmatrix(:,1) == day,:);
            epochs = runmatrix(:,2);
            
            daycells = animcells(animcells(:,1) == day,[2,3]);
            
            for c = 1:size(daycells,1)
                flag = 0;
                cell = daycells(c,:);
                
                
                si = cellinfo{day}{epochs(1)}{cell(1)}{cell(2)}.SI;
                rightfr = []; leftfr = [];
                lfields = []; rfields = [];
                for ep = 1:length(epochs)
                    epoch = epochs(ep);
                    
                    if length(linfields{day}{epoch}{cell(1)}) >= cell(2)
                        lfield = linfields{day}{epoch}{cell(1)}{cell(2)}{1}(:,5);
                        lfields = stack(lfields,lfield');
                        
                        rfield = linfields{day}{epoch}{cell(1)}{cell(2)}{3}(:,5);
                        rfields = stack(rfields,rfield');
                        
                    else
                        continue
                    end
                    
                end
                
                lfield = nanmean(lfields,1);
                rfield = nanmean(rfields,1);
                
                maxlength = min(length(lfield), length(rfield));
                fields = [lfield(1:maxlength),rfield(1:maxlength)];
                [~,peakfield] = max(fields);
                
                %determine if peak is on right or left
                if peakfield > maxlength
                    peakfield = peakfield-maxlength;
                    field = rfield;
                else
                    field = lfield;
                end
                
                
                %find window around peak
                if peakfield == 1
                    win(1) = 1;
                    win(2) = find(field(peakfield:maxlength) <= thresh, 1, 'first') + peakfield -1 ;
                elseif peakfield == maxlength
                    win(2) = peakfield;
                    win(1) = find(field(1:peakfield) <= thresh, 1, 'last');
                else
                    if any(field(1:peakfield) <= thresh)
                        win(1) = find(field(1:peakfield) <= thresh, 1, 'last');
                    else
                        win(1) = 1;
                    end
                    
                    if any(field(peakfield:end) <= thresh)
                        win(2) = find(field(peakfield:end) <= thresh, 1, 'first') + peakfield -1 ;
                    else
                        win(2) = maxlength;
                    end
                end
                
                
                lfr = lfield(win(1):win(2));
                rfr = rfield(win(1):win(2));
                
                
                %end
                mnl = nanmean(lfr);
                mnr = nanmean(rfr);
                
                fieldsi = (mnl - mnr)/(mnl + mnr);
                %                 if abs(fieldsi) == 1
                %                     keyboard
                %                 end
                
                fieldcomps = [fieldcomps; si, fieldsi];
            end
        end
    end
    
    [bad,~] = find(isnan(fieldcomps));
    fieldcomps(bad,:) = [];
    
    figure,
    plot(fieldcomps(:,1), fieldcomps(:,2), 'k.','MarkerSize',20);
    hold on
    fit = polyfit(fieldcomps(:,1), fieldcomps(:,2),1);
    plot([-1, 1], polyval(fit,[min(fieldcomps(:,1)), max(fieldcomps(:,1))]))
    axis([-1 1 -1 1]);
    xlabel('Odor Selectivity Index');
    ylabel('Spatial Trajectory Selectivity Index');
    
    [CC,p] = corrcoef(fieldcomps(:,1), fieldcomps(:,2));
    R = CC(1,2)
    p = p(1,2)
    
    text(0.5,0.7,{['R = ' num2str(R)],;['p = ' num2str(p)]})
    title(region)
    
    figfile = [figDir, 'PlaceFields\SI-PlaceFieldCorrelation_', region];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);
end