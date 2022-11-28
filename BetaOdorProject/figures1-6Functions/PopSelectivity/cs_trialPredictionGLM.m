%cs_predictChoice_GLM
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%regions = {'CA1','PFC'};
load([topDir,'AnalysesAcrossAnimals\npCells_CA1']);
regions = {'CA1','PFC'};
c1 = npCells;
load([topDir,'AnalysesAcrossAnimals\npCells_PFC'])
c2 = npCells;
allcells = [c1;c2];
                    
predict = 'Choice';

celltype = 'all';
switch celltype
    case 'all'
        celltag = '';
    case 'selective'
        celltag = '_selective';
end

binsize = 0.1;
maxwin = 1;

wins = binsize:binsize:maxwin;

for r = 1:length(regions)
    region = regions{r};
    
    Acc_mn = {};
    Acc = {};
    Acc_shuff = {};
    for w = 1:length(wins)
        window = wins(w);
        
        predAcc = [];
        predAcc_shuff = [];
        for a = 1:length(animals)
            animal = animals{a};
            topDir = cs_setPaths();
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            runeps = cs_getRunEpochs(animDir, animal, 'odorplace');
            days = unique(runeps(:,1));
            for day = days'
                disp(['Doing ',num2str(window),'s window prediction, animal ',animal,' day ',num2str(day)])
                if strcmp(region,'all')
                    cells = allcells(ismember(allcells(:,[1 2]),[a,day],'rows'),:);
                    
                    [fract_correct, fract_correct_shuff] = cs_predictChoice_GLM_v2(animal, day, window,cells);
                else
                    switch predict
                        case 'Choice'
                            [fract_correct, fract_correct_shuff] = cs_predictChoice_GLM(animal, day, region, window,celltype);
                        case 'Outcome'
                            [fract_correct, fract_correct_shuff] = cs_predictOutcome_GLM(animal, day, region, window,celltype);
                    end
                end
                
                
                predAcc = [predAcc;fract_correct];
                predAcc_shuff = [predAcc_shuff; fract_correct_shuff'];
            end
        end
        Acc_mn{end+1} = mean(predAcc);
        Acc{end+1} = predAcc;
        Acc_shuff{end+1} = predAcc_shuff;
        
        save([topDir, 'AnalysesAcrossAnimals\predict',predict,'GLM_',region,celltag],'Acc_mn','Acc','Acc_shuff','wins');
        
    end
end