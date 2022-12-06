function f = setfiltereegtetrodepairs(f, filterStringArray, options)
% f = setfiltereegtetrodepairs(f, filterStringArray)
% For each epoch in the filter F, this function finds pairs of tetrodes that
% that satisfy each of the filter conditions in the 1st and 2nd
% element of the filterStringArray.  
%
% The final list is a double array with one row for each pair and four columns:
% dataset epoch tet1 tet2.
%
% Note that this list is put into f.eegdata
%
% The syntax for each filterString is defined in EVALUATEFILTER.m. 
% The animal and desired epochs for the filter need to be predefined. 
% Assumes that each animal's data  folder contains a file ending in 
% 'tetinfo.mat' that contains a cell structure with information about each cell.

if (length(filterStringArray) ~= 2) 
    error('filterStringArray must have 2 elements');
end

for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the cells'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering the cells'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    tetinfo = loaddatastruct(datadir,animalprefix,'tetinfo');

    for i = 1:length(f(an).epochs)
        if isempty(f(an).epochs{i})
            f(an).data{i} = [];
        end
        for j = 1:size(f(an).epochs{i},1)
	    for k = 1:length(filterStringArray)
		tmplist{k} = evaluatefilter(tetinfo{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}, filterStringArray{k});
		tmpind{k} = 1:size(tmplist{k}, 1);
	    end
	    % combine the lists to get all combinations
	    % filterStringArray
	    [i1 i2] = meshgrid(tmpind{1}, tmpind{2});
	    celllist = [tmplist{1}(i1,:) tmplist{2}(i2,:)];
	    % get rid of identical elements
	    valid = (sum(celllist(:,1) == celllist(:,2), 2) ~= 2);
	    celllist = celllist(valid,:);
	    if (length(celllist))
		% check for pairs listed twice because the order was reversed
		tmpcelllist = [celllist(:,2) celllist(:,1)];
		[c i1 i2] = intersect(celllist, tmpcelllist, 'rows');
		valid = ones(size(celllist,1), 1);
		% for each member of tmpcelllist that is in celllist, we
		% need to delete the correponding member of celllist if
		% it's index is larger than the index in tmpcelllist
		valid(i2) = (i2 < i1);
		f(an).eegdata{i}{j} = celllist(logical(valid),:);
	    else
		f(an).eegdata{i}{j} = [];
	    end
        end
    end
end
