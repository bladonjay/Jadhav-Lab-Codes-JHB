function out = loaddatastruct(animaldir, animalprefix, datatype, days)

% out = loaddatastruct(animaldir, animalprefix, datatype)
% out = loaddatastruct(animaldir, animalprefix, datatype, days)
%
% Load the components of a data cell array of combines then into one
% variable.  Datatype is a string with the base name of the files.  If DAYS
% is omitted, all files are loaded.  Otherwise, only the files for the
% specified days will be included.
%
% Ryan -  To warn the user that a clear request for data returned void, I'm
% adding a warning; hope this will help future users save some time
% troubleshooting code. It's one of the most common problems in execution
% I've noticed, rooted in calls to this function with an incorrect prefix
% path. Not always immediately apparent when it happens, especially when
% nested in the parent function. See line 25 for change.

if (nargin < 4)
    days = [];
end
out = [];
% if strcmp(datatype, 'lowgamma') == 1
%     animaldir = [animaldir, 'EEG\'];
% end
datafiles = dir([animaldir,animalprefix, datatype,'*']);

% Ryan -- code to warn user if no files found
if isempty(datafiles)
	errormessage = ['Looked in %s for animal %s''s ' ...
				'datatype ''''%s'''' and failed to find it' '\n\n'];
	warning(errormessage,animaldir,animalprefix,datatype);
end

for i = 1:length(datafiles)
    if isempty(days)
        load([animaldir,datafiles(i).name]);
        eval(['out = datavaradd(out,',datatype,');']);
    else
        s = datafiles(i).name;
        fileday = str2num(s(strfind(s,datatype)+length(datatype):strfind(s,'.')-1));  %get the experiment day from the filename
        if (isempty(fileday))|(ismember(fileday,days))
            load([animaldir,datafiles(i).name]);
            if strcmp(datatype,'ripplesep1'), % file ripplesep1 is already loaded
                datatype = 'ripples';
            end
            eval(['out = datavaradd(out,',datatype,');']);
        end
    end
end

%--------------------------------------
function out = datavaradd(origvar, addcell)

out = origvar;
for i = 1:length(addcell)
    if (~isempty(addcell{i}))
        out{i} = addcell{i};
    end
end
        
