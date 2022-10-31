% convert indicies and data to matricies
% find all unique rows in matrix based on columns 1, 3, 4
%for each unique row, find original index rows that share those numbers in
%1 3 4
% for those rows in data matrix, average or concatenate and enter into new
% data structure, one row for each unique cell with epochs combined (new
% index - [day tet cell] only) 
% 
function [gathered] = cs_gatherepochs(datadir, prefix, datatype, fieldpath, fieldstocombine, newfieldnames, fieldstokeep)

alldata = load([datadir,prefix,datatype,'.mat']);

alldata = alldata.(datatype).(fieldpath);

rawindices = vertcat(alldata.index);
daysonly = unique(rawindices(:,[1,3,4]),'rows');

for i = 1:length(fieldstocombine)
    fieldname = fieldstocombine{i};
    
    %rawdata = vertcat(alldata.(fieldname));

    for j = 1:length(daysonly)
        samedaytetcell = find((rawindices(:,1) == daysonly(j,1)) & (rawindices(:,3) == daysonly(j,2)) & (rawindices(:,4) == daysonly(j,3)));
    
        combined = [];
        for t = 1:length(samedaytetcell)
             epochindex = samedaytetcell(t);
             combined = [combined; alldata(epochindex).(fieldname)];
        end
             
%         datatocombine = rawdata(samedaytetcell,:);
%         combinedepochsdata = vertcat(datatocombine); %maybe concatenate instead? 
%         
        newfieldname = newfieldnames{i};
    
        combineddata(j).(newfieldname) = combined;
        combineddata(j).index = daysonly(j,:);
    
    
        for f = 1:length(fieldstokeep)
            combineddata(j).(fieldstokeep{f}) = alldata(1).(fieldstokeep{f});
        
        end
    end
    
    gathered = combineddata;
    
end


        
    
    
    
    