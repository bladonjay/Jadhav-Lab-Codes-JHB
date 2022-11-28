
animals={'CS31','CS33','CS34'}; %CHANGE DATADIR PATH IN 'animaldef.mat' FOR DIFFERENT COMPUTERS
    
topDir = 'F:\Data\OdorPlaceAssociation\'; %home computer
dataDir = [topDir,'AnalysesAcrossAnimals\'];

    runepochfilter = 'isequal($environment, ''odorplace'')';
    cellfilter = '((strcmp($area, ''CA1'') && ($numspikes > 100))  ||  (strcmp($area, ''PFC'') && ($numspikes > 100) ) )';
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6} }; 
    riptetfilter = '(isequal($descrip, ''riptet''))';
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= 0) & (abs($linearvel) >= 5))', 3, 'headdir',1}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    %timefilter = {{'DFTFsj_getvelpos', '(($absvel >= 5))'}};
    %tetfilter = '(strcmp($area, ''CA1''))';
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
     psf = createfilter('animal',animals,'epochs',runepochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator', iterator);
    % psf = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'iterator', iterator);
    
    % Set analysis function
    % ---------------------
    pmf_all = setfilterfunction(psf, 'DFAcssj_openfieldrate_leftrightall', {'spikes', 'linpos', 'pos', 'runTrialBounds'}, 'appendindex' ,1 , 'binsize', 1, 'std', 2, 'trialtype', 2);    
%     pmf_left = setfilterfunction(psf, 'DFAcs_linposfield', {'spikes', 'linpos', 'runTrialBounds'}, 'appendindex' ,1 , 'binsize', 1, 'std', 2, 'trialtype', 1);    
%     pmf_right = setfilterfunction(psf, 'DFAcs_linposfield', {'spikes', 'linpos', 'runTrialBounds'}, 'appendindex' ,1 , 'binsize', 1, 'std', 2, 'trialtype', 0); 
    % Run analysis 
    % ------------
    pmf_all = runfilter(pmf_all);  % Place Field Map
    %pmf_left = runfilter(pmf_left);
    %pmf_right = runfilter(pmf_right);
    
    % tetrode metadata
%  hpc   = [2 3 4 17 18 19 20 27 28 29 30 31 32];
%  hpcRef=  1;
%  pfc   = [7 8 9 11 12 13 14 15 16 21 22 23 25 26];
%  pfcRef=   10;
%  ob = [24];
 
 totalcells = 0;
 for a = 1:length(pmf_all)
     animal = placeFieldData(a).animal{1, 1};
     data = placeFieldData(a).output{1, 1};
     
     cellindex = vertcat(data.index);
     cellindexnoepochs = cellindex(:,[1,3,4]);
     uniquecells = unique((cellindexnoepochs),'rows');
     
     for c = 1:length(uniquecells)
            cell = uniquecells(c,:);
            cellstocombine = data(ismember(cellindexnoepochs,cell, 'rows'));
            
                for i = 1:length(cellstocombine)
                    sizes(i,:) = size(cellstocombine(i).smoothedspikerate);
                end
                
                height = min(sizes(:,1));
                width = min(sizes(:,2));
                
                clear sizes
                
                for i = 1:length(cellstocombine)
                    cellstocombine(i).smoothedspikerate = cellstocombine(i).smoothedspikerate(1:height,1:width);
                    cellstocombine(i).smoothedoccupancy = cellstocombine(i).smoothedoccupancy(1:height,1:width);
                end
                
                newdata(totalcells+c).animal = animal;
                newdata(totalcells+c).index = cell;
                newdata(totalcells+c).smoothedspikerate = mean(cat(3,cellstocombine.smoothedspikerate),3);
                newdata(totalcells+c).smoothedoccupancy = mean(cat(3,cellstocombine.smoothedoccupancy),3);
     
     
     end
     
     totalcells = length(newdata);
   
     
 end
 
 cellregions = {'CA1','PFC'};
 for r = 1:length(cellregions)
     region = cellregions{r};
     
     totalcells = 0;
     for a = 1:length(animals)
         animal = animals{a};
         animalstrs = {newdata.animal};
         animalcells = newdata(strcmp(animalstrs,animal));


         load([topDir,animal,'Expt\',animal,'_direct\',animal,'tetinfo.mat'])
         tetfilter = ['((strcmp($area, ''',region,''')))'];
         tets = evaluatefilter(tetinfo{1,1}{1,2}, tetfilter);
         
         cellinds = vertcat(animalcells.index);
         celltets = cellinds(:,2);
         good = (ismember(celltets,tets, 'rows'));
         cellsinregion(totalcells+1:(totalcells+sum(good))) = animalcells(good);
         totalcells = length(cellsinregion);
     end
     
     eval(['placeFields_',region,'= cellsinregion;']);
     clear cellsinregion
 end
 
 save([dataDir, 'placeFieldData.mat'],'placeFields_CA1','placeFields_PFC');
         
     
 