 regions = {'CA1','PFC','OB','TC'};
 animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%regions = {'PFC','OB'};
%animals = {'CS44'};
[topDir] = cs_setPaths();

runfilter = 'strcmp($environment,''odorplace'') || strcmp($environment,''novelodor'') || strcmp($environment,''novelodor2'')';

do_wrtgnd = 1;
fpass = [0 15];
movingwin = [1000 20]/1000;

for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    load([animDir, animal,'tetinfo.mat']);
    
    task = loaddatastruct(animDir, animal,'task');
    ep = evaluatefilter(task, runfilter);

    days = unique(ep(:,1));
    epochs = unique(ep(:,2));
    
    for r = 1:length(regions)
        region = regions{r};
        disp(['Doing ', animal,' ', region]);
        if strcmp(region,'CA1') || strcmp(region,'PFC')
        tetfilter = ['strcmp($area,''',region,''') && $numcells>0'];
        else
            tetfilter = ['strcmp($area,''',region,''')'];
        end
        alltets = evaluatefilter(tetinfo,tetfilter);
        
        tets = unique(alltets(:,3));
        
        sjcs_baselinespecgram(animal, days, epochs, tets, do_wrtgnd,'fpass',fpass,'movingwin',movingwin);
    end
end