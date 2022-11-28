%CS 5/7/2020

%calculates PCs for left and right trials. Plot to visualize difference between population

%cs_PCA(region, win, binsize, selectiveonly)
%close all
%cs_PCA('CA1',[0 1], 0.05, 1)

%% --- Params
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
[topDir, figDir] = cs_setPaths();
figdir = [figDir,'PCA\'];

selective = 1;
allCells = [];
for r = 1:length(regions)
    region = regions{r};
    
    if selective == 1
        selstr = '_selective';
        load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region]);
        allCells= [allCells;selectivecells];
    else
        selstr = '_np';
        load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
        allCells= [allCells;npCells];
    end
end

%% --- Create Column Vectors
disp('creating column vectors')

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    nosepokeFR = loaddatastruct(animDir, animal,'nosepokeFR');
    
    cells = allCells(allCells(:,1) == a,[2:4]);
    
    days = unique(cells(:,1));
    
    for day = days'
        
        
        daycells = cells(cells(:,1) == day,:);
        
        if size(daycells) >=3
        runEps = cs_getRunEpochs(animDir, animal,'odorplace',day);
        epochs = unique(runEps(:,2));
        %Gather trig times
        left_master = [];
        right_master = [];
        for c = 1:size(daycells)
            cell = daycells(c,:);
            
            leftfr_cell = [];
            rightfr_cell = [];
            for ep = epochs'
                
                [cl, cr, il, ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                
                leftfr = nosepokeFR{day}{ep}{cell(2)}{cell(3)}([cl]);
                rightfr = nosepokeFR{day}{ep}{cell(2)}{cell(3)}([cr]);
                
                leftfr_cell = [leftfr_cell,leftfr'];
                rightfr_cell = [rightfr_cell,rightfr'];
            end
            
            left_master = [left_master; leftfr_cell];
            right_master = [right_master; rightfr_cell];
        end
        
        LRcombined = [left_master, right_master];
        disp(['calculating pca ',animal,' ',num2str(day)])
        [coeff,score,~,~,explained] = pca(LRcombined');
        new = score(:,1:2)*coeff(1:2,:);
        newL = new(1:size(left_master,2),:);
        newR = new(size(left_master,2)+1:end,:);
        

        figure,
        plot(newL(:,1),newL(:,2),'b.')
        hold on
        plot(newR(:,1),newR(:,2),'r.')
        xlabel('PC1');
        ylabel('PC2');
        
        
        close
        end
    end
end


%%
