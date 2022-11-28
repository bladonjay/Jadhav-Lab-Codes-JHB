
function cs_getNPerrors(animals)

[topDir,~] = cs_setPaths();

for a = 1:length(animals)
    animal = animals{a};
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    
    cd(dataDir)
    
    epochs = cs_getRunEpochs(dataDir, animal,'odorplace');
%    
    days = unique(epochs(:,1));
    
    for d = 1:length(days)
        day = days(d);
        daystring = getTwoDigitNumber(day);
        
        
        runEps = epochs(epochs(:,1) == day, 2);
        
        load([animal,'DIO', daystring, '.mat']);
        numEpochs = length(runEps);
        
        for e = 1:length(runEps)
            ep = runEps(e);
            
            
            %%  Get DIO times and states for reward wells, odor solenoids, nosepoke, and buzzer
            
            Din1times = dio{1,day}{1,ep}{1,1}.time(logical(dio{1,day}{1,ep}{1,1}.state));
            Din1hits = ones([length(Din1times),1]);
            leftwellnans = nan([length(Din1hits),2]);
            LeftWell = [Din1times, leftwellnans, Din1hits];
            clear('Din1times','Din1hits','leftwellnans');
            
            Din2times = dio{1,day}{1,ep}{1,2}.time(logical(dio{1,day}{1,ep}{1,2}.state));
            Din2hits = ones([length(Din2times),1])+1;
            rightwellnans = nan([length(Din2hits),2]);
            RightWell = [Din2times, rightwellnans, Din2hits];
            clear('Din2times','Din2hits','rightwellnans');
            
            Reward = sortrows([LeftWell; RightWell]);
            

            %%%% Buzzer is now on Dout5 - use for all animals after CS31 %%%%
            
            if strcmp(animal, 'CS31') || strcmp(animal, 'CS42')
                Dout9times = dio{1,day}{1,ep}{1,25}.time(logical(dio{1,day}{1,ep}{1,25}.state));
                Dout9hits = ones([length(Dout9times),1]);
                buzzernans = nan([length(Dout9hits),1]);
                Buzzer = [Dout9times, buzzernans, Dout9hits, buzzernans];
                clear('Dout9times','Dout9hits','buzzernans');
                
            else
                
                Dout5times = dio{1,day}{1,ep}{1,21}.time(logical(dio{1,day}{1,ep}{1,21}.state));
                Dout5hits = ones([length(Dout5times),1]);
                buzzernans = nan([length(Dout5hits),1]);
                Buzzer = [Dout5times, buzzernans, Dout5hits, buzzernans];
                clear('Dout5times','Dout5hits','buzzernans');
                
            end
            
            Din5times = dio{1,day}{1,ep}{1,5}.time;
            Din5states = double(dio{1,day}{1,ep}{1,5}.state);
            npNans = nan([length(Din5states),2]);
            NP = [Din5times, Din5states, npNans];
            clear('Din5times','Din5hits','NPnans');
            
            
            %%
            %create matrix with np state, which odor was dispensed, when buzzer
            %went off, and which reward well was trigered.
            
            attempts = sortrows([NP; Buzzer; Reward]);
            
            %for each np, determine whether next hit was buzzer (good
            %trial) or reward (np error). If np error, save the timestamp.
            
            
            nperrs = [];
            npInds = find((attempts(:,2)==1));
            for n = 1:length(npInds)
                
                npind = npInds(n);
                %what was the next trigger?
                if n ~= length(npInds)
                nextnp = npInds(n+1);
                else
                    nextnp = size(attempts,1)+1;
                end
                nextbuzzer = find(~isnan(attempts(npind+1:end,3)),1,'first')+npind;
                nextreward = find(~isnan(attempts(npind+1:end,4)),1,'first')+npind;
                
                if nextreward < nextbuzzer & nextreward < nextnp
                    nperrs = [nperrs; attempts(n,1)];
                end
                
                npErrors{day}{ep} = nperrs;
            end
            
        end
        save([dataDir,animal, 'npErrors', daystring],'npErrors');
        clear npErrors
        
    end
    
end
