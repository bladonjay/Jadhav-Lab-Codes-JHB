function cs_coherencePos_nponly(animals, topDir, regions, figDir, approach)

%approach = 1, times before nosepoke
%approach = 0, times after nosepoke
%approach = 2, window around nosepoke
close all
for r = 1:length(regions)
    region = regions{r};
    
    allData = [];
    figure, hold on;
    for a = 1:length(animals)
        animal = animals{a};
        
        coherencefiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'coherence',region,'*']);
        posfiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'pos*']);
        odorfiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'odorTriggers*']);
        
        for d = 1:length(coherencefiles)
            day = d; 
            load([topDir,animal,'Expt\',animal,'_direct\',coherencefiles(d).name])
            load([topDir,animal,'Expt\',animal,'_direct\',posfiles(d).name])
            load([topDir,animal,'Expt\',animal,'_direct\',odorfiles(d).name])
            
            coherence = coherence{1,d};
            pos = pos{1,d};
            odorTriggers = odorTriggers{1,d};
            epochs = find(~cellfun(@isempty, coherence));
            
            for ep = 1:length(epochs)
                
                epoch = epochs(ep);
                betaFreqs = find((coherence{1,epoch}.freqs > 15) & (coherence{1,epoch}.freqs < 30));
                meanBetaCoh = mean(coherence{1,epoch}.coherence(betaFreqs,:))';
                
                
                
                cohTimes = coherence{1,epoch}.time';
                posTimes = pos{1,epoch}.data(:,1);
                newPos = round(pos{1,epoch}.data(:,2:3)); %round to 1cm bins
                
                %frame rate for pos is much slower than samp rate for coherence.
                %Interpolate to get new beta power. 
                interpCoh = interp1(cohTimes, meanBetaCoh, posTimes);
                
                cohPos = [posTimes, newPos, interpCoh];
                
                trigs = odorTriggers{1,epoch}.allTriggers;
                goodtimes = [];
                
                    for t = 1:length(trigs)
                        if approach == 1
                            goodtimes = [goodtimes; posTimes((posTimes < trigs(t)) & (posTimes > (trigs(t) -0.5)))];
                        elseif approach == 0
                            goodtimes = [goodtimes; posTimes((posTimes < (trigs(t) +1.5)) & (posTimes > trigs(t)))];
                        elseif approach == 2
                            times = posTimes((posTimes < (trigs(t) +1.5)) & (posTimes > trigs(t) - 1));
                            posplot = cohPos(ismember(cohPos(:,1),times),:);
                            plot(posplot(:,3),posplot(:,2),'k');
                            goodtimes = [goodtimes; times];
                            
                        end
                    end
                
                cohPos = cohPos(ismember(cohPos(:,1),goodtimes),:);

                meanCohPos = accumarray(cohPos(:,2:3),cohPos(:,4),[],@mean); 
                
                
                [xsize] = max([size(meanCohPos,1), size(allData,1)]);
                [ysize] = max([size(meanCohPos,2), size(allData,2)]);
                
                %put new data in the middle of existing matrix, padded with
                %zeros on edges.
                meanCoh_adjusted = zeros(xsize,ysize);
                xstart = round((xsize - size(meanCohPos,1))/2)+1;   
                ystart = round((ysize - size(meanCohPos,2))/2)+1;
                meanCoh_adjusted(xstart:(xstart+size(meanCohPos,1)-1), ystart:(ystart+size(meanCohPos,2)-1)) = meanCohPos;
                
                allData_adjusted = zeros(xsize,ysize,size(allData,3));
                xstart = round((xsize - size(allData,1))/2)+1;   
                ystart = round((ysize - size(allData,2))/2)+1;
                allData_adjusted(xstart:(xstart+size(allData,1)-1), ystart:(ystart+size(allData,2)-1),:) = allData;
                
                
                allData_adjusted(:,:,size(allData_adjusted,3)+1) = meanCoh_adjusted;
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
    imagesc(flipud(smoothedData))
    colormap(jet)
    title(region)
    colorbar
    
    
    if approach == 1
        descriptstring = 'approach';
    elseif approach == 0
        descriptstring = 'exit';
    elseif approach == 2
        descriptstring = 'window';
    end
    
    figfile = [figDir,'1f_BetaCohPos_',region,'_',descriptstring];

    print('-djpeg', figfile);
    print('-dpdf', figfile);
    saveas(gcf,figfile,'fig');
end
