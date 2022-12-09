%Binomial test for whether animal turns away from nosepoke in same
%direction as reward

%load animal pos, np window files, task files, and odortriggers (to determine which
%direction animal ended up going)

[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
timebins = 40;
binomial_all = [];
ci = 1; %Do correct/incorrect separately
binomial_c = [];
binomial_i = [];
savefig=0;

prcnt_all = [];
for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir, animal, 'Expt\',animal, '_direct\'];
    
    dayeps = cs_getRunEpochs(animDir, animal,'odorplace');
    nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow');
    pos = loaddatastruct(animDir, animal,'pos');
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    task = loaddatastruct(animDir,animal,'task');
    
    days = unique(dayeps(:,1));
    
    %files = dir([animDir, animal, 'nosepokeWindow*']);
    
    for day = days'

        daystr = getTwoDigitNumber(day);
        
%         load([animDir, animal, 'pos',daystr,'.mat']);
%         load([animDir, animal, 'odorTriggers',daystr,'.mat'])
%         load([animDir, animal, 'task',daystr,'.mat']);
        if strcmpi(animal,'CS41') && day<3
            epochs = 1; % claire collapsed epochs for this animal only
            %dio{1,day}{1}=dio{1,day}; % have to add back in a cell for epoch
        else
            epochs = dayeps(dayeps(:,1) == day,2);
        end

        for epoch = epochs'
           
            [c_left, c_right, i_left, i_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
            
            
            if ci == 0
                %% --- All Trials Together ---% 
                %in this case, "left trials" are those in which the animal went to
                %the left reward well, regardless of which odor was presented. The
                %time is nosepoke EXIT time (to start looking at pos)
                leftinds = sort([c_left; i_right]);
                rightinds = sort([i_left; c_right]);
                
                lefttrials = nosepokeWindow{day}{epoch}(leftinds,2);
                righttrials = nosepokeWindow{day}{epoch}(rightinds,2);
                if ~isfield(task{day}{epoch},'linearcoord')
                    continue
                else
                    npPos = task{day}{epoch}.linearcoord{1}(1,1,1); %only need xpos
                    
                    epochpos = pos{day}{epoch}.data(:,[1:3]);
                    
                    %diffs = (epochpos(2:end,1)-epochpos(1:end-1,1));
                    
                    lefttimeinds = lookup(lefttrials, epochpos(:,1));
                    righttimeinds = lookup(righttrials, epochpos(:,1));
                    
                    %remove trials too close to end of epoch
                    lefttimeinds = lefttimeinds(lefttimeinds+40 < size(epochpos,1));
                    righttimeinds = righttimeinds(righttimeinds+40 < size(epochpos,1));
                    
                    
                    for L = 1:length(lefttimeinds)
                        ind = lefttimeinds(L);
                        xvals = epochpos(ind:ind+timebins,2);
                        diffs = xvals- npPos;
                        %find furthest point from NP
                        [~,maxind] = max(abs(diffs));
                        point = diffs(maxind);
                        
                        %for left turns, max distance will be negative, and for
                        %right turns, max distance will be positive
                        if point < 0
                            binomial(end+1) = 1;
                        else
                            binomial(end+1) = 0;
                        end
                    end
                    
                    
                    for R = 1:length(righttimeinds)
                        ind = righttimeinds(R);
                        xvals = epochpos(ind:ind+timebins,2);
                        diffs = xvals- npPos;
                        %find furthest point from NP
                        [~,maxind] = max(abs(diffs));
                        point = diffs(maxind);
                        %for left turns, max distance wills be negative, and for
                        %right turns, max distance will be positive
                        if point < 0
                            binomial(end+1) = 0;
                        else
                            binomial(end+1) = 1;
                        end
                    end
                    
                end
                
                binomial_all = [binomial_all,binomial];
                prcnt = sum(binomial)/length(binomial);
                prcnt_all = [prcnt_all;prcnt];
            
            %% --- Correct/Incorrect Separate --- %%    
            else
                %in this case, "left trials" are those in which the animal went to
                %the left reward well, regardless of which odor was presented. The
                %time is nosepoke EXIT time (to start looking at pos)
%                 leftinds = sort([c_left; i_right]);
%                 rightinds = sort([i_left; c_right]);
                
                %CORRECT
                lefttrials = nosepokeWindow{day}{epoch}(c_left,2);
                righttrials = nosepokeWindow{day}{epoch}(c_right,2);
                if ~isfield(task{day}{epoch},'linearcoord')
                    continue
                else
                    npPos = task{day}{epoch}.linearcoord{1}(1,1,1); %only need xpos
                    
                    epochpos = pos{day}{epoch}.data(:,[1:3]);
                    
                    %diffs = (epochpos(2:end,1)-epochpos(1:end-1,1));
                    
                    lefttimeinds = lookup(lefttrials, epochpos(:,1));
                    righttimeinds = lookup(righttrials, epochpos(:,1));
                    
                    %remove trials too close to end of epoch
                    lefttimeinds = lefttimeinds(lefttimeinds+40 < size(epochpos,1));
                    righttimeinds = righttimeinds(righttimeinds+40 < size(epochpos,1));
                    
                    
                    for L = 1:length(lefttimeinds)
                        ind = lefttimeinds(L);
                        xvals = epochpos(ind:ind+timebins,2);
                        diffs = xvals- npPos;
                        %find furthest point from NP
                        [~,maxind] = max(abs(diffs));
                        point = diffs(maxind);
                        
                        %for left turns, max distance wills be negative, and for
                        %right turns, max distance will be positive
                        if point < 0
                            binomial_c(end+1) = 1;
                        else
                            binomial_c(end+1) = 0;
                        end
                    end
                    
                    
                    for R = 1:length(righttimeinds)
                        ind = righttimeinds(R);
                        xvals = epochpos(ind:ind+timebins,2);
                        diffs = xvals- npPos;
                        %find furthest point from NP
                        [~,maxind] = max(abs(diffs));
                        point = diffs(maxind);
                        %for left turns, max distance wills be negative, and for
                        %right turns, max distance will be positive
                        if point < 0
                            binomial_c(end+1) = 0;
                        else
                            binomial_c(end+1) = 1;
                        end
                    end
                    
                end
                
                %INCORRECT
                righttrials = nosepokeWindow{day}{epoch}(i_left,2); %incorrect left means he ultimately went to the right
                lefttrials = nosepokeWindow{day}{epoch}(i_right,2);
                if ~isfield(task{day}{epoch},'linearcoord')
                    continue
                else
                    npPos = task{day}{epoch}.linearcoord{1}(1,1,1); %only need xpos
                    
                    epochpos = pos{day}{epoch}.data(:,[1:3]);
                    
                    %diffs = (epochpos(2:end,1)-epochpos(1:end-1,1));
                    
                    lefttimeinds = lookup(lefttrials, epochpos(:,1));
                    righttimeinds = lookup(righttrials, epochpos(:,1));
                    
                    %remove trials too close to end of epoch
                    lefttimeinds = lefttimeinds(lefttimeinds+40 < size(epochpos,1));
                    righttimeinds = righttimeinds(righttimeinds+40 < size(epochpos,1));
                    
                    
                    for L = 1:length(lefttimeinds)
                        ind = lefttimeinds(L);
                        xvals = epochpos(ind:ind+timebins,2);
                        diffs = xvals- npPos;
                        %find furthest point from NP
                        [~,maxind] = max(abs(diffs));
                        point = diffs(maxind);
                        
                        %for left turns, max distance wills be negative, and for
                        %right turns, max distance will be positive
                        if point < 0
                            binomial_i(end+1) = 1;
                        else
                            binomial_i(end+1) = 0;
                        end
                    end
                    
                    
                    for R = 1:length(righttimeinds)
                        ind = righttimeinds(R);
                        xvals = epochpos(ind:ind+timebins,2);
                        diffs = xvals- npPos;
                        %find furthest point from NP
                        [~,maxind] = max(abs(diffs));
                        point = diffs(maxind);
                        %for left turns, max distance wills be negative, and for
                        %right turns, max distance will be positive
                        if point < 0
                            binomial_i(end+1) = 0;
                        else
                            binomial_i(end+1) = 1;
                        end
                    end
                    
                end
                
            end
        end
    end
end
disp('Done');

if ci == 0
    
    mn = mean(prcnt_all);
    sem = stderr(prcnt_all);
prob = binofit((length(binomial) - sum(binomial)),length(binomial),0.5);
pie([prob, 1-prob]);

else
    
save([topDir, 'AnalysesAcrossAnimals\turnDirection'],'binomial_c','binomial_i');
    prob_c = binofit((length(binomial_c) - sum(binomial_c)),length(binomial_c),0.5);
subplot(1,2,1)
pie([prob_c, 1-prob_c]);
title('Correct Trials');

y = binopdf(sum(binomial_c),length(binomial_c),0.5)

prob_i = binofit((length(binomial_i) - sum(binomial_i)),length(binomial_i),0.5);
subplot(1,2,2)
pie([prob_i, 1-prob_i]);
title('Incorrect Trials')

figfile = [figDir, 'Behavior\turnDirectionFromNP'];
           if savefig==1 
    print('-djpeg', figfile);
    print('-dpdf', figfile);
           end
y = binopdf(sum(binomial_i),length(binomial_i),0.5)

% [p,z] = ztestprop(sum(binomial_c), length(binomial_c)-sum(binomial_c), length(binomial_c), length(binomial_c))
% [p,z] = ztestprop(sum(binomial_i), length(binomial_i)-sum(binomial_i), length(binomial_i), length(binomial_i))

end