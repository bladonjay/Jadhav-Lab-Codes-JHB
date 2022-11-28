%cs_getNPerrors
clear

topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39'};

for a = 1:length(animals)
    animal = animals{a};
    
    dataDir =([topDir,animal,'Expt\',animal,'_direct\']);
    epochs = cs_getRunEpochs(dataDir, animal, 'odorplace');
    days = unique(epochs(:,1));
    
    for d = 1:length(days)
        day = days(d);
        daystring = getTwoDigitNumber(day);
        
        load([animal,'DIO', daystring, '.mat']);
        
        
        runEps = epochs(find(epochs(:,1) == day),2);
        
        for e = 1:length(runEps)
            errors = [];
            epoch = runEps(e);
            
            %         Din5times = dio{1,day}{1,epoch}{1,5}.time;
            %         Din5states = double(dio{1,day}{1,epoch}{1,5}.state);
            %         allnps = Din5times(logical(Din5states));
            
            Din1times = dio{1,day}{1,epoch}{1,1}.time(logical(dio{1,day}{1,epoch}{1,1}.state));
            Din1hits = ones([length(Din1times),1]);
            leftwellnans = nan([length(Din1hits),2]);
            LeftWell = [Din1times, leftwellnans, Din1hits];
            clear('Din1times','Din1hits','leftwellnans');
            
            Din2times = dio{1,day}{1,epoch}{1,2}.time(logical(dio{1,day}{1,epoch}{1,2}.state));
            Din2hits = ones([length(Din2times),1])+1;
            rightwellnans = nan([length(Din2hits),2]);
            RightWell = [Din2times, rightwellnans, Din2hits];
            clear('Din2times','Din2hits','rightwellnans');
            
            Reward = sortrows([LeftWell; RightWell]);
            
            if strcmp(animal, 'CS31')
                Dout9times = dio{1,day}{1,epoch}{1,25}.time(logical(dio{1,day}{1,epoch}{1,25}.state));
                Dout9hits = ones([length(Dout9times),1]);
                buzzernans = nan([length(Dout9hits),1]);
                Buzzer = [Dout9times, buzzernans, Dout9hits, buzzernans];
                clear('Dout9times','Dout9hits','buzzernans');
                
            else
                
                Dout5times = dio{1,day}{1,epoch}{1,21}.time(logical(dio{1,day}{1,epoch}{1,21}.state));
                Dout5hits = ones([length(Dout5times),1]);
                buzzernans = nan([length(Dout5hits),1]);
                Buzzer = [Dout5times, buzzernans, Dout5hits, buzzernans];
                clear('Dout5times','Dout5hits','buzzernans');
                
            end
            
            Din5times = dio{1,day}{1,epoch}{1,5}.time(logical(dio{1,day}{1,epoch}{1,5}.state));
            Din5hits = ones([length(Din5times),1]);
            npNans = nan([length(Din5hits),2]);
            NP = [Din5times, Din5hits, npNans];
            clear('Din5times','Din5hits','NPnans');
            
            attempts = sortrows([NP; Buzzer; Reward]);
            
            npInds = find(attempts(:,2)==1);
            for i = 1:length(npInds)
                ind = npInds(i);
                
                nextnp = find(~isnan(attempts(ind+1:end,2)),1,'first')+ind;
                nextbuzzer = find(~isnan(attempts(ind+1:end,3)),1,'first')+ind;
                prevbuzzer = find(~isnan(attempts(1:ind,3)),1,'last');
                prevreward = find(~isnan(attempts(1:ind,4)),1,'last');
                nextreward = find(~isnan(attempts(ind+1:end,4)),1,'first')+ind;
                
                if ~isempty(nextreward)
                    if (isempty(nextnp) || nextreward < nextnp) && (isempty(nextbuzzer) || nextbuzzer > nextreward)...
                            && (isempty(prevbuzzer) || isempty(prevreward) || prevbuzzer < prevreward)
                       
                        trigger = attempts(ind,1);
                        errors = [errors; trigger];
                        
                    end
                end
                
            end
            
            %%%
            
            nperrors{day}{epoch} = errors;
        end
        
        save([dataDir,animal,'nperrors',daystring,'.mat'],'nperrors');
        clear nperrors
    end
end