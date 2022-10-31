%find tetrode with largest difference between baseline and after-trigger
%beta
function betatet = cs_findBetaTet_v2(animal, day, region)
topDir = cs_setPaths();
daystr = getTwoDigitNumber(day);
animDir = [topDir,animal,'Expt\',animal,'_direct\'];

dayepochmatrix = cs_getRunEpochs(animDir, animal, 'odorplace',day);
runeps = unique(dayepochmatrix(:,2));
nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow',day);
tetinfo = loaddatastruct(animDir, animal,'tetinfo');

allwins = cellfun(@(x) x,nosepokeWindow{day},'UniformOutput',0);
allwins = vertcat(allwins{:});

tetfilter = ['(isequal($area, ''',region,'''))'];
tets = evaluatefilter(tetinfo{day}, tetfilter);
tets = unique(tets(:,2));


zscores = zeros(length(tets),1);
for t = 1:length(tets)
    tet = tets(t);
    tetstr = getTwoDigitNumber(tet);
    
    betamag = [];
    betatime = [];
    for e = 1:length(runeps)
        epoch = runeps(e);
        epstr = getTwoDigitNumber(epoch);
        %                
        try
            load([animDir,'EEG\',animal,'beta',daystr,'-',epstr,'-',tetstr,'.mat'])
        catch
            continue
        end
        
        betamag = [betamag;double(beta{day}{epoch}{tet}.data(:,3))];
        time = geteegtimes(beta{day}{epoch}{tet})';
        betatime = [betatime;time];
    end
    betaaftertrig = mean(betamag(isExcluded(betatime,allwins)));
    pretrigs = [allwins(:,1) - (allwins(:,2)-allwins(:,1)), allwins(:,1)];
    betapretrig = betamag(isExcluded(betatime,pretrigs));
    baseline = mean(betapretrig);
    baselinestdev = std(betapretrig);
    
    zscores(t) = (betaaftertrig - baseline) ./ baselinestdev;
end



[val,ind] = max(zscores);
betatet = tets(ind);

% if no tets had increase in beta, use
% tet with most cells
 if val <= 0
     betatet = cs_getMostCellsTet(animal,day,epoch,region);
 end

end