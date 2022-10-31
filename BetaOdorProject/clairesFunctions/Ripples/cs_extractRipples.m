%cs_extractRipples
%script that loops over all animals and extracts ripples from eeg
%uses riptets (defined in cs_riptettag, all CA1 tets with cells on them)
min_suprathresh_duration = 0.015;
nstd = 3;

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

topDir = cs_setPaths();

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir,animal,'Expt',filesep,animal,'_direct',filesep];
    tetinfo = loaddatastruct(animDir, animal,'tetinfo');
    days = cs_getRunEpochs(animDir, animal,'odorplace');
    days = unique(days(:,1));
    for d = days'
        tetfilter = 'strcmp($descrip,''riptet'')';
        riptets = evaluatefilter(tetinfo{d},tetfilter);
        riptets = unique(riptets(:,2));
        disp(['Doing ', animal, ' day ',num2str(d)])
        sj_extractripples(animDir, animal, d, riptets, min_suprathresh_duration, nstd);
    end
end

cs_rippletimes;