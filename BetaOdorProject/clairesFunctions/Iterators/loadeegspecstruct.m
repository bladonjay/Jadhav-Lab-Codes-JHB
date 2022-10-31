function out = loadeegspecstruct(animaldir, animalprefix, datatype, days, tet)

% out = loadeegstruct(animaldir, animalprefix, datatype, days, epochs, tetrodes)
%
% Load the components of an eeg cell array of combines then into one
% variable. 
%	datatype is a string with the base name of the files (e.g. 'theta') 
% 	days specifies the list of days to load 
% 	epochs specifies the list of epochs to load 
% 	tetrodes specifies the list of tetrodes to load 
%
%	if epochs is empty, all epochs (1-10) are loaded
%	if tetrodes is empty, all tetrodes (1-50) are loaded
%
% Be aware that eeg files tend to be large, so loading a large number of them
% is likely to lead to problems

if isempty(tet)
    tet = 1:50;
end


out = [];
% create the list of data files 

for d = days
        for t = tet
            if isunix
                fname = sprintf('%s/EEGSpec/%s%s%02d-%02d.mat', animaldir,...
                animalprefix, datatype, d, t);
            else
                fname = sprintf('%s\\EEGSpec\\%s%s%02d-%02d.mat', animaldir,...
                animalprefix, datatype, d, t);
            end
            try
                load(fname);
                filename = datatype(1:end-3);
                    if strcmp(filename, 'eeggndspecfl')
                        filename = 'eeggndspec';
                    end
                eval(['out = datavaradd(out,',filename,');']);
            catch
                disp('Error in loadeegstruct');
                keyboard
            end
        end
end


%--------------------------------------
function out = datavaradd(origvar, addcell)

out = origvar;
for d = 1:length(addcell)
    for e = 1:length(addcell{d})
	for t = 1:length(addcell{d}{e})
	    if (~isempty(addcell{d}{e}{t}))
		out{d}{e}{t} = addcell{d}{e}{t};
	    end
	end
    end
end
