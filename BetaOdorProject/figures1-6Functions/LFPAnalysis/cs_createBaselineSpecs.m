function cs_createBaselineSpecs(animals, regions, do_wrtgnd, fpass, movingwin, tapers, trialtype)
%this script loops over animals and creates baseline spectrograms based on
%specified parameters. Eliminates having to run each animal separately.
% this runs the function over ALL TETRODES

% example run:
%cs_createBaselineSpecs({'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'}, {'CA1'}, 1, [0 15], [2000 25]/1000, [2 1])

%cs_createBaselineSpecs({'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'}, {'CA1'}, 0, [0 300], [100 10]/1000, [2 1])

%animals = {'CS41','CS42'};
%regions = {'OB','TC'};

%do_wrtgnd = 1;

topDir = cs_setPaths();


for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
    runepochs = cs_getRunEpochs(animDir, animal, trialtype);
    days = unique(runepochs(:,1));
    
    for r = 1:length(regions)
        region = regions{r};
        %tetfilter = ['isequal($area,''',region,''') & strcmp($descrip,''riptet'')'];
        tetfilter = ['isequal($area,''',region,''')& ($numcells > 1)'];
        for day = days'
            tets = evaluatefilter(tetinfo{day},tetfilter);
            tets = unique(tets(:,2));
            
            epochs = (runepochs(runepochs(:,1) == day,2))';
            sjcs_baselinespecgram(animal, day, epochs, tets, do_wrtgnd,'fpass',fpass,'movingwin',movingwin,'tapers',tapers);
        end
    end
end