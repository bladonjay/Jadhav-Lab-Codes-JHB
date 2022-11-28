function cs_betaPowerPos_nponly(animals, topDir, regions, figDir, approach)

%approach = 1 for approach, approach = 0 for exit 

for r = 1:length(regions)
    region = regions{r};
    
    allData = [];
    
    for a = 1:length(animals)
        animal = animals{a};
        
        betafiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'meanBeta',region,'*']);
        posfiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'pos*']);
        odorfiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'odorTriggers*']);
        
        for d = 1:length(betafiles)
            day = d; 
            load([topDir,animal,'Expt\',animal,'_direct\',betafiles(d).name])
            load([topDir,animal,'Expt\',animal,'_direct\',posfiles(d).name])
            load([topDir,animal,'Expt\',animal,'_direct\',odorfiles(d).name])
            
            meanBeta = meanBeta{1,d};
            pos = pos{1,d};
            odorTriggers = odorTriggers{1,d};
            epochs = find(~cellfun(@isempty, meanBeta));
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                betaTimes = meanBeta{1,epoch}(:,1);
                posTimes = pos{1,epoch}.data(:,1);
                newPos = round(pos{1,epoch}.data(:,2:3));
                
                %frame rate for pos is much slower than samp rate for LFP.
                %Interpolate to get new beta power. 
                interpBeta = interp1(betaTimes, meanBeta{1,epoch}(:,2), posTimes);
                
                betaPos = [posTimes, newPos, interpBeta];
                
                trigs = odorTriggers{1,epoch}.allTriggers;
                goodtimes = [];
                    for t = 1:length(trigs)
                        if approach == 1
                            goodtimes = [goodtimes; posTimes((posTimes < trigs(t)) & (posTimes > (trigs(t) -1)))];
                        else
                            goodtimes = [goodtimes; posTimes((posTimes < (trigs(t) +1.5)) & (posTimes > trigs(t)))];
                        end
                    end
                betaPos = betaPos(ismember(betaPos(:,1),goodtimes),:);
                %for each position bin, find the average beta power
                
                meanPower = accumarray(betaPos(:,2:3),betaPos(:,4),[],@mean); 
                
                [xsize] = max([size(meanPower,1), size(allData,1)]);
                [ysize] = max([size(meanPower,2), size(allData,2)]);
                
                %put new data in the middle of existing matrix, padded with
                %zeros on edges.
                meanPower_adjusted = zeros(xsize,ysize);
                xstart = round((xsize - size(meanPower,1))/2)+1;   
                ystart = round((ysize - size(meanPower,2))/2)+1;
                meanPower_adjusted(xstart:(xstart+size(meanPower,1)-1), ystart:(ystart+size(meanPower,2)-1)) = meanPower;
                
                allData_adjusted = zeros(xsize,ysize,size(allData,3));
                xstart = round((xsize - size(allData,1))/2)+1;   
                ystart = round((ysize - size(allData,2))/2)+1;
                allData_adjusted(xstart:(xstart+size(allData,1)-1), ystart:(ystart+size(allData,2)-1),:) = allData;
                
                
                allData_adjusted(:,:,size(allData_adjusted,3)+1) = meanPower_adjusted;
                allData = allData_adjusted;
                
                
                        
            end
                
        end
    end
    
    
    allData = mean(allData,3);
    std = 2; g = gaussian2(std,(4*std)); 
    smoothedData = filter2(g,(allData)); % 
    
    [goodx, goody] = find(smoothedData);
    smoothedData = smoothedData(min(goodx):max(goodx),min(goody):max(goody));
    smoothedData = padarray(smoothedData, [8 8], 0, 'both');
    
    figure,
    imagesc(smoothedData)
    colormap(jet)
    colorbar
    title(region)
    
    if approach == 1
        descriptstring = 'approach';
    else
        descriptstring = 'exit';
    end
    
    figfile = [figDir,'1f_BetaPowerPos_',region,'_',descriptstring];

    print('-djpeg', figfile);
    print('-dpng', figfile);
    saveas(gcf,figfile,'fig');
end

%for each region, for each animal, for each epoch(?) load pos, load betapower

% interpolate (?) pos times with beta times

%plot in a similar way to place fields - occupancy normalized? 