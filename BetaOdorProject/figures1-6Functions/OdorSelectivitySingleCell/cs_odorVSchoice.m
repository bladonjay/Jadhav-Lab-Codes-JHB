%cs_odorVSchoice
%separates sells that exhibit selectivity into those that are purely odor
%selective (i.e. don't care about incorrect/correct) and those that are
%choice selective

regions = {'CA1','PFC'};
topDir = cs_setPaths;


for r = 1:length(regions)
    odorcells = [];
    choicecells = [];
    region = regions{r};
    load([topDir,'AnalysesAcrossAnimals\selectivityData_',region]);
    si_correct = selectivityData.SI_correct;
    si_incorrect = selectivityData.SI_incorrect;
    for c = 1:size(si_correct,1)
        if sign(si_correct(c)) == sign(si_incorrect(c))
            odorcells = [odorcells; selectivityData.cellinds(c,:)];
        else 
            choicecells = [choicecells; selectivityData.cellinds(c,:)];
        end
    end
save([topDir,'AnalysesAcrossAnimals\odorCells_',region],'odorcells');
save([topDir,'AnalysesAcrossAnimals\choiceCells_',region],'choicecells');

end

