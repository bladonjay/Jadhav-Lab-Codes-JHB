%1/20/20
%Plot rasters for selective cells
%loads data from cs_selectivityRaster


%Good cells for paper:
% 3_3_30_1 - best
% 7_4_2_9 
% 5_1_4_1 
%PFC:
% 8_1_15_2
% 7_5_9_1 
% 2_2_22_2  - best

%% Params
clear
[topDir, figDir] = cs_setPaths;
close all
regions = {'CA1','PFC'};

savefig = 1;

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

leftcolor = rgb('BrightRoyalBlue');
rightcolor = rgb('Crimson');
cellstoplot = [3 3 30 1;2 2 22 2]; %top row = CA1, bottom row = PFC. Leave empty to plot all cells
%cellstoplot = [];
%% Get data and plot
for r = 1:length(regions)
    region = regions{r};
    
    load([topDir, 'AnalysesAcrossAnimals\rasterData_',region]);
    
    inds = rasterData.inds;
    
    if ~isempty(cellstoplot)
        cellind = find(ismember(rasterData.inds,cellstoplot(r,:),'rows'));
        rasterData.inds = rasterData.inds(cellind,:);
        rasterData.data{1} = rasterData.data{cellind};
    end
    
    for c = 1:size(rasterData.inds,1)
        index = rasterData.inds(c,:);
        animal = animals{index(1)};
        
        figure, set(gcf,'Position',[500 300 1200 400]);
        subplot(1,2,1)
        
        %correct trials
        leftspikes = rasterData.data{c}.cleft;
        leftspikes = leftspikes(randperm(length(leftspikes)));
        
        rightspikes =  rasterData.data{c}.cright;
        
       rightspikes = rightspikes(randperm(length(rightspikes)));   
          
        Ly = [];  
        for L = 1:length(leftspikes)
            if isempty(leftspikes{L})
                continue
            else
            Ly = [Ly, repmat(L,1,length(leftspikes{L}))];
            end
        end
        notempty = find(~cellfun(@isempty, leftspikes));
        leftspikes = {leftspikes{notempty}};
        leftspikes = cell2mat(leftspikes');
        
%         
%         
        Ry = [];
        for R = 1:length(rightspikes)
            if isempty(rightspikes{R})
                continue
            else
            Ry = [Ry, repmat(R+max(Ly),1,length(rightspikes{R}))];
            end
        end
        notempty = find(~cellfun(@isempty, rightspikes));
        rightspikes = {rightspikes{notempty}};
        rightspikes = cell2mat(rightspikes');
        
        scatter(leftspikes, Ly, 100, '.', 'MarkerEdgeColor', leftcolor);
        hold on
        scatter(rightspikes, Ry,100,'.','MarkerEdgeColor',rightcolor);
        
        plot([0 0],[0 max(Ry)],'k:')
        axis([-0.4 1 0 max(Ry)])

        Ry_c = Ry;
        xlabel('Time from odor onset (s)');
        ylabel('Trial Number');

        
        %incorrect trials
        subplot(1,2,2)

        leftspikes = rasterData.data{c}.ileft;
        rightspikes =  rasterData.data{c}.iright;
        
        if isempty(leftspikes) || isempty(rightspikes)
            continue
        end
          
        Ly = [];  
        for L = 1:length(leftspikes)
            if isempty(leftspikes{L})
                continue
            else
            Ly = [Ly, repmat(L,1,length(leftspikes{L}))];
            end
        end
         notempty = find(~cellfun(@isempty, leftspikes));
         leftspikes = {leftspikes{notempty}};

        leftspikes = cell2mat(leftspikes');
        
        
        Ry = [];
        for R = 1:length(rightspikes)
            if isempty(rightspikes{R})
                continue
            else
            Ry = [Ry, repmat(R+max(Ly),1,length(rightspikes{R}))];
            end
        end
        notempty = find(~cellfun(@isempty, rightspikes));
        rightspikes = {rightspikes{notempty}};

        rightspikes = cell2mat(rightspikes');
        
        scatter(leftspikes, Ly, 100, '.', 'MarkerEdgeColor', leftcolor);
        hold on
        scatter(rightspikes, Ry,100,'.','MarkerEdgeColor',rightcolor);
        
        plot([0 0],[0 max(Ry_c)],'k:')
        axis([-0.4 1 0 max(Ry_c)])
        xlabel('Time from odor onset (s)');
        ylabel('Trial Number');
        
        figtitle = ['Raster_',region,'_',animal,'_',num2str(index(2)),'_',num2str(index(3)),'_',num2str(index(4))];
        figfile = [figDir,'OdorSelectivity\',figtitle];
        
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        
    end
    close all
end