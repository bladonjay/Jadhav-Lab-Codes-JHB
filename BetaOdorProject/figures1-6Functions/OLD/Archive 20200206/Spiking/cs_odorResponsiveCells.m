%cs_odorResponsiveCells
%include all cells that have a CHANGE in FR during odor presentation. Cells
%can be task responsive if they are inhibited, and even if they are not
%selective. Spiking during NP does not necessarily mean cells are
%task-responsive if there is no specific change

%compares between window before NP to window after NP.

clear 
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
%regions = {'CA1'};
baselineLength = 0.5;

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
    region = regions{r};
    
    odorRespCells =[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['isequal($area,''',region,''')'];
%         cellfilter = ['isequal($type, ''pyr'')']
%         cellfilter = ['isequal($area,''',region,''')']
        cells = evaluatefilter(cellinfo,cellfilter);

        cells = unique(cells(:,[1 3 4]),'rows');
        
        days = unique(cells(:,1));
        for d = 2:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'nosepokeWindow',daystr,'.mat'])
            load([animDir,animal,'pos',daystr,'.mat'])
            runeps = find(~cellfun(@isempty,nosepokeWindow{day}));
            
           
            for c = 1:size(daycells,1)

                npfr = [];
                ctrlfr = [];
                totnptime = 0;
                totaltrigs = 0;
                psthspikes = [];
                cell = daycells(c,:);
                npposall = [];
                posall = [];
                for ep = 1:length(runeps)
                    epoch = runeps(ep);

                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                            epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);

                            trigs = nosepokeWindow{day}{epoch};
                            %baselinetrigs = [trigs(:,1)- (trigs(:,2)-trigs(:,1)), trigs(:,1)];
                            baselinetrigs = [trigs(:,1)- baselineLength, trigs(:,1)];

                            %totnptime = totnptime + sum(trigs(:,2)-trigs(:,1));
                            
                            totaltrigs = totaltrigs + size(trigs,1);
                            for tr = 1:size(trigs,1)
                                winlength = trigs(tr,2) - trigs(tr,1);
                            winspikesraw = epspikes(isExcluded(epspikes, trigs(tr,:)));
                            npfr = [npfr; length(winspikesraw)/winlength];
                            winspikes = winspikesraw - trigs(tr,1);

                            ctrlspikesraw = epspikes(isExcluded(epspikes,baselinetrigs(tr,:)));
                            ctrlfr = [ctrlfr; length(ctrlspikesraw)/winlength];
                            ctrlspikes = ctrlspikesraw- trigs(tr,1);
                            psthspikes = [psthspikes; winspikes; ctrlspikes];
                            
                            %also get pos, check where spikes are actually
                            %occuring
                            epspikes = [winspikesraw;ctrlspikesraw];
                            npposind = lookup(epspikes,pos{cell(1)}{epoch}.data(:,1));
                            nppos = pos{cell(1)}{epoch}.data(npposind,[2,3]);
                            npposall = [npposall;nppos];
                            posall = [posall; pos{cell(1)}{epoch}.data(:,[2,3])];
                            end
                        
                    end

                end
                if isempty(npfr) || (~any(npfr) && ~any(ctrlfr))
                    continue
                end
                [~,test] = ttest2(npfr, ctrlfr)
                if test <= 0.05
                    odorRespCells = [odorRespCells; a, cell];
                end
                counts = histcounts(psthspikes,[-baselineLength:0.05:1.5+0.05]);
                psth = (counts/totaltrigs) / 0.05;
                psth = smoothdata(psth,'gaussian',5);
                figure,
                plot([-baselineLength:0.05:1.5],psth);
                hold on
                plot([0 0],[0 max(psth)],'r--');
                text(-1, 1, ['p = ',num2str(test)]);
                
                figure, plot(posall(:,1),posall(:,2));
                hold on
                plot(npposall(:,1),npposall(:,2),'r.');
                
            end
        end
    end
    
    save([dataDir,'odorResponsiveCells_',region,'.mat'], 'odorRespCells')
end