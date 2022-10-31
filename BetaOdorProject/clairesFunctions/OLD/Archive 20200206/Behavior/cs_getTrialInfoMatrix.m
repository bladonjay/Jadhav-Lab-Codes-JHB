%make matrix with trial info 
%[odor time, l/r (0/1), correct/incorrect (0/1)];

function trialinfo = cs_getTrialInfoMatrix(animal, day, epoch)

topDir = cs_setPaths();

animDir = [topDir, animal,'Expt\',animal,'_direct\'];
odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
trigstruct = odorTriggers{day}{epoch};
[correct_left, correct_right, ~, incorrect_right] = cs_getSpecificTrialTypeInds(trigstruct);

times = trigstruct.allTriggers;
lr = zeros(length(times),1);
lr([correct_right;incorrect_right]) = 1;
lr = lr + 6;

ci = zeros(length(times),1);
ci([correct_left;correct_right]) = 1;

trialinfo = [times, lr, ci];
