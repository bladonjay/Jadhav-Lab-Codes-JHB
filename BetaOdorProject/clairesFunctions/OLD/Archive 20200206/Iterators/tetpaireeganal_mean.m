function f = tetpaireeganal_mean(f)
% f = singlecelleeganal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.
% Also the function call appends any options in the f().function.options{}
% cell array.
%
% This iterator should be used when you want to use perform an analysis where
% you select a single eeg channel and a set of single cells and wish to
% perform an analysis that combines the eeg and single cell data.
%
% Note that for any specified eeg variables, only the data for the tetrode
% specified by the eegfilter will be loaded and passed to the function.
%
% Each function call is for one cell (or one tetrode), and it is assumed that
% the function's first input is the index to the cell or tetrode ([day epoch tetrode
% cell]).  The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

eegvar = [];

%iterate through all animals
for an = 1:length(f)
    %find all unique epochs to analyze for the current animal
    animaldir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    disp(['Doing ',animalprefix]);
    totalepochs = [];
    for g = 1:length(f(an).epochs)
        totalepochs = [totalepochs; f(an).epochs{g}];
    end
    totaldays = unique(totalepochs(:,1)); %get all of the days across groups
    
    %load all the variables that the function requires
    loadstring = [];
    for i = 1:length(f(an).function.loadvariables)
        if (~iseegvar(f(an).function.loadvariables{i}))
            eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
        end
        loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        
        for e = 1:size(f(an).epochs{g},1)
                 
            for c = 1:size(f(an).eegdata{g}{e},1)
                tet1index = [f(an).epochs{g}(e,:) f(an).eegdata{g}{e}(c,1)];
                tet2index = [f(an).epochs{g}(e,:) f(an).eegdata{g}{e}(c,2)];
                % load the eeg data for this cell
                eegload = {};
                for i = 1:length(f(an).function.loadvariables)
                    if (iseegvar(f(an).function.loadvariables{i}))
                        eval([f(an).function.loadvariables{1},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tet1index(1), tet1index(2), tet1index(3));']);
                        eegload{tet1index(1)}{tet1index(2)}{tet1index(3)} = eeg{tet1index(1)}{tet1index(2)}{tet1index(3)};
                        eval([f(an).function.loadvariables{1},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tet1index(1), tet1index(2), tet2index(3));']);
                        eegload{tet2index(1)}{tet2index(2)}{tet2index(3)} = eeg{tet2index(1)}{tet2index(2)}{tet2index(3)};
                    else 
                       eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tet1index(1));']);
                    end
                end
                
                excludeperiods = f(an).excludetime{g}{e};
                
                eeg = eegload;
                
                %run the designated function: fout = fname(sindex, eindex, var1, var2, ..., option1, option2, ...)
                eval(['fout = ',f(an).function.name,'(tet1index, tet2index, excludeperiods,', loadstring, 'foptions{:});']);
                
                %save the function output in the filter variable.  Allows numeric or struct outputs
                if isstruct(fout)
                    if (isempty(f(an).output) | (length(f(an).output) < g))
                        f(an).output{g}(1) = fout;
                    else
                        f(an).output{g}(end+1) = fout;
                    end
                elseif isnumeric(fout)
                    if ((isempty(f(an).output)) | (length(f(an).output) < g))
                        f(an).output{g} = [];
                    end
                    if (size(fout,2) >= 1) %Changed by MCarr to allow N by 1 numeric functions (previously had been inconsistent)
                        f(an).output{g} = [f(an).output{g}; fout];
                    else
                        error(['In calling ', f(an).function.name, ': Numeric function outputs must be 1 by N.  Use a structure output for more complicated outputs']);
                    end
                else
                    error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
                end
            end
            
            data = cat(3, f(an).output{g}.epochCoherence); 
            meanData = mean(data,3);
            f(an).output{g}(e).meanData = meanData;  %take the mean over all tetrode pairs for the epoch. Do this here rather than at the end to save space. 
            f(an).output{g}(e).epoch = f(an).epochs{g}(e,:);
            f(an).output{g}(e).params.binlength = f(an).output{1, 1}.binlength;
            f(an).output{g}(e).params.freqs = f(an).output{1, 1}.freqs;
            f(an).output{g}(e).params.time = fout.time ;
            [f(an).output{g}.epochCoherence] = deal([]);
        end
    end
    f(an).output{g} = rmfield(f(an).output{g},'epochCoherence');
    f(an).output{g} = rmfield(f(an).output{g},'binlength'); 
    f(an).output{g} = rmfield(f(an).output{g},'index'); 
    f(an).output{g} = rmfield(f(an).output{g},'freqs'); 
    f(an).output{g} = rmfield(f(an).output{g},'time'); 
end




