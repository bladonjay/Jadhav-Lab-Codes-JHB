
%Plot nice PSTH and rasters for specified odor selective cells.
%enter cells to plot as a matrix. 

% stored in "selectiveCells_(region)" files in AnalysesAcrossAnimals

%% Params
[topDir, figDir] = cs_setPaths;
close all
region = 'CA1';

savefig = 1;

load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region,'.mat']);
cellstoplot = selectivecells;

win = [0.5 1];
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
stdev = 4;
g = gaussian(stdev,(2*stdev));
leftcolor = rgb('BrightRoyalBlue');
rightcolor = rgb('Crimson');

%%

for c = 1:size(cellstoplot,1)
    cindex = cellstoplot(c,:);
    animal = animals{cindex(1)};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    
    day = cindex(2);
    daystr = getTwoDigitNumber(day);
    
     load([animDir,animal,'spikes',daystr,'.mat'])
     load([animDir,animal,'odorTriggers',daystr,'.mat'])
     load([animDir,animal,'cellinfo.mat'])
     
     %region = cellinfo{day}{2}{cindex(3)}{cindex(4)}.area;
     
     runeps = find(~cellfun(@isempty,odorTriggers{day}));
     
     cf = '$numspikes > 0';
     allcells = evaluatefilter(cellinfo,cf);
     
     %epochs where that cell actually spiked
     celleps = allcells(ismember(allcells(:,[1 3 4]), cindex(2:4),'rows'),2);
     
     runeps = celleps(ismember(celleps,runeps));
     
  
     alltrigs = []; lefttrigs = []; righttrigs = [];
     for ep = 1:length(runeps)
              epoch = runeps(ep);
              alltrigs = [alltrigs; odorTriggers{day}{epoch}.allTriggers];
              [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
              lefttrigs = [lefttrigs; alltrigs(correct_left)];
              righttrigs = [righttrigs; alltrigs(correct_right)];
     end
     
    
     allspikes = [];
      for ep = 1:length(runeps)
          epoch = runeps(ep);

          if ~isempty(spikes{cindex(2)}{epoch}{cindex(3)}{cindex(4)})
                allspikes = [allspikes; spikes{cindex(2)}{epoch}{cindex(3)}{cindex(4)}.data(:,1)];
          end


      end
      
       leftspikes = []; 
       rightspikes = [];
     
       for t = 1:length(lefttrigs)
            trig = lefttrigs(t);
                trigwin = [trig-win(1), trig+win(2)];
                winspikes = allspikes(allspikes > trigwin(1) & allspikes <= trigwin(2));
                winspikes = winspikes- trig;

                leftspikes{t} = winspikes';     
       end
       
%        tmp = cell2mat(leftspikes);
%        leftspikecount = histcounts(tmp,bins);
%        avgspikecount = leftspikecount./length(lefttrigs); %avg num spikes per trial
%        lfr = avgspikecount./binsize; %convert to hz
       
       %leftPSTH = filter(g,1,lfr); 
       %plot(leftPSTH);
       %hold on
       %plot(lfr)
       %leftPSTH = smoothdata(lfr,'gaussian',5);
       % plot(leftPSTH);

       %lfr = [lfr(1),lfr,lfr(end)];
       %leftPSTH = filter(g,1,lfr); 
       %plot(leftPSTH(2:end-1));

       for t = 1:length(righttrigs)
            trig = righttrigs(t);
                trigwin = [trig-win(1), trig+win(2)];
                winspikes = allspikes(allspikes > trigwin(1) & allspikes <= trigwin(2));
                winspikes = winspikes- trig;
                

                rightspikes{t} = winspikes';     
       end
        
%         tmp = cell2mat(rightspikes);
%         rightspikecount = histcounts(tmp, bins);
%         avgspikecount = rightspikecount./length(righttrigs);
%         rfr = avgspikecount./binsize; %convert to hz
%         %rightPSTH = filter2(g,rfr); 
%         rightPSTH = smoothdata(rfr,'gaussian',5);
        %plot(leftPSTH);
        
        %figure, set(gcf,'Position',[800 40 700 600]);
        
        %leftspikes = cell2mat(leftspikes);
        
         
         notempty = find(~cellfun(@isempty, leftspikes));
         leftspikes = {leftspikes{:,notempty}};
          
        Ly = [];  
        for L = 1:length(leftspikes)
            Ly = [Ly, repmat(L,1,length(leftspikes{L}))];
        end
        leftspikes = cell2mat(leftspikes);
        
        notempty = find(~cellfun(@isempty, rightspikes));
        rightspikes = {rightspikes{:,notempty}};
        
        Ry = [];
        for R = 1:length(rightspikes)
            Ry = [Ry, repmat(R,1,length(rightspikes{R}))];
        end
        rightspikes = cell2mat(rightspikes);
        
        subplot(2,1,1)
        scatter(leftspikes, Ly, 100, '.', 'MarkerEdgeColor', leftcolor);
        axis([-win(1) win(2) 0 L+1])
        xlabel('Time from Odor Onset (seconds)')
        ylabel('Trials')
        
        subplot(2,1,2)
        scatter(rightspikes, Ry,100,'.','MarkerEdgeColor',rightcolor);
        axis([-win(1) win(2) 0 R+1])
        xlabel('Time from Odor Onset (seconds)')
        ylabel('Trials')
        
        
        %% Plot
%         hold on
%         plot(bins(1:end-1), leftPSTH,  'LineWidth',3.5, 'Color', leftcolor);
%         plot(bins(1:end-1), rightPSTH, 'LineWidth',3.5, 'Color', rightcolor);
% 
%         maxfr = max([leftPSTH,rightPSTH]);
% 
%         plot([0 0], [0 maxfr+1], 'k--','LineWidth',2);
%         ylabel('Firing Rate (Hz)')
%         xlabel('Time from Odor Onset (seconds)')
%         set(gca,'fontsize',25);
%         axis([-win(1) win(2) 0 maxfr+1])
       
    
    
    
        if savefig == 1
            figtitle = ['Raster_',region,'_',animal,'_',num2str(cindex(2)),'_',num2str(cindex(3)),'_',num2str(cindex(4))];
            figfile = [figDir,'OdorSelectivity\',figtitle];
            
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            
            saveas(gcf,figfile,'fig');
        end

        %close all
  
end