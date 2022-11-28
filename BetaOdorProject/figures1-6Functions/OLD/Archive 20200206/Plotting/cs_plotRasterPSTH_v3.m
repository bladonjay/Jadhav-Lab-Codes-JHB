% cs_plotRasterPSTH_v3

win = [0.5 1];
binsize = 0.05;
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
bins = (-win(1):binsize:win(2));

stdev = 3;
g = gaussian(stdev,(3*stdev));

for a = 1:length(animals)
    animal = animals{a};
    
    load([topDir,animal,'Expt\',animal,'_direct\',animal,'cellinfo.mat']);
    
    cellfilter = ['(strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective''))'];
    cells = evaluatefilter(cellinfo,cellfilter);
    
    days = unique(cells(:,1));
    
    for d = 1:length(days)
        day = days(d); daystr = getTwoDigitNumber(day);
        
        daycells = cells(cells(:,1) == day,:);
        daycells = unique(daycells(:,[1 3 4]),'rows'); %remove epochs
        
        
        load([topDir,animal,'Expt\',animal,'_direct\',animal,'spikes',daystr,'.mat'])
        load([topDir,animal,'Expt\',animal,'_direct\',animal,'odorTriggers',daystr,'.mat'])
        
         runeps = find(~cellfun(@isempty,odorTriggers{day}));
         
         
         alltrigs = []; lefttrigs = []; righttrigs = [];
         for ep = 1:length(runeps)
                  epoch = runeps(ep);
                  alltrigs = [alltrigs; odorTriggers{day}{epoch}.allTriggers];
                  [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                  lefttrigs = [lefttrigs; alltrigs(correct_left)];
                  righttrigs = [righttrigs; alltrigs(correct_right)];
         end
         
         for c = 1:size(daycells,1)
              
              cell = daycells(c,:);
              
              allspikes = [];
              for ep = 1:length(runeps)
                  epoch = runeps(ep);
                  
                  if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        allspikes = [allspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];
                  end
                  
                  
              end
              
              
               leftspikes = []; 
               rightspikes = [];
               
               for t = 1:length(lefttrigs)
                    trig = lefttrigs(t);
                        trigwin = [trig-win(1), trig+win(2)];
                        winspikes = allspikes(allspikes > trigwin(1) & allspikes <= trigwin(2));
                        winspikes = winspikes- trig;
                        %bins = (trig-win(1):binsize:trig+win(2));
                        %binspikecount = histcounts(winspikes,bins);
                        
                        leftspikes = [leftspikes; winspikes];     
               end
               
               leftspikecount = histcounts(leftspikes,bins);
               avgspikecount = leftspikecount./length(lefttrigs);
               lfr = avgspikecount./binsize; %convert to hz
               leftPSTH = filter2(g,lfr); 
               
               for t = 1:length(righttrigs)
                    trig = righttrigs(t);
                        trigwin = [trig-win(1), trig+win(2)];
                        winspikes = allspikes(allspikes > trigwin(1) & allspikes <= trigwin(2));
                        winspikes = winspikes- trig;
                        %bins = (trig-win(1):binsize:trig+win(2));
                        %binspikecount = histcounts(winspikes,bins);
                        
                        rightspikes = [rightspikes; winspikes];     
                end
                
                rightspikecount = histcounts(rightspikes, bins);
                avgspikecount = rightspikecount./length(righttrigs);
                rfr = avgspikecount./binsize; %convert to hz
                rightPSTH = filter2(g,rfr); 
                
                
                figure, 
                set(gcf,'Position',[1800 250 560 420]);
                %suptitle([animal,' Cell ', num2str(cell)])
                subplot(3, 1, 1), hold on
                for i = 1:size(leftspikes,1)
                    plot(leftspikes(i,:), i, 'g.')
                end

                plot([0 0], [0 (length(leftspikes))], 'k--');
                axis( [-win(1) win(2) 0 length(leftspikes)])

                subplot(3, 1, 2), hold on
                for i = 1:size(rightspikes,1)
                    plot(rightspikes(i,:), i ,'m.')
                end
                plot([0 0], [0 (length(rightspikes))], 'k--');
                axis( [-win(1) win(2) 0 length(rightspikes)])

                subplot(3, 1, 3)
                hold on
                plot(bins(1:end-1), leftPSTH, 'g');
                plot(bins(1:end-1), rightPSTH, 'm');
                
                
                tmp = [leftPSTH, rightPSTH];
                ymax = max(tmp);
                
                axis( [-win(1) win(2) 0 ymax+1])
                plot([0 0], [0 ymax+1], 'k--');
        
        
        

       pause;
                
                
         end
    end
end


