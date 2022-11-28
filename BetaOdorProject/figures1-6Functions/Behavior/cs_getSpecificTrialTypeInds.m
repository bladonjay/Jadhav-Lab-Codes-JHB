function [correct_left, correct_right, incorrect_left, incorrect_right,prelearn,postlearn] = cs_getSpecificTrialTypeInds(trigstruct)
%load trigstruct from odorTriggers from desired day and epoch (do not load
%raw odorTriggers file)

%outputs trial index for each trial type

[~,correct_left] = intersect(trigstruct.allTriggers, intersect(trigstruct.correctTriggers, trigstruct.leftTriggers));
[~,incorrect_left] = intersect(trigstruct.allTriggers, intersect(trigstruct.incorrectTriggers, trigstruct.leftTriggers));
[~,correct_right] = intersect(trigstruct.allTriggers, intersect(trigstruct.correctTriggers, trigstruct.rightTriggers));
[~,incorrect_right] = intersect(trigstruct.allTriggers, intersect(trigstruct.incorrectTriggers, trigstruct.rightTriggers));
if isfield(trigstruct,'prelearn')
    [~,prelearn] = intersect(trigstruct.allTriggers, trigstruct.prelearn);
    [~,postlearn] = intersect(trigstruct.allTriggers, trigstruct.postlearn);
else
    prelearn = [];
    postlearn = [];
end

check = sum([length(correct_left), length(incorrect_left), length(correct_right), length(incorrect_right)]);

if check ~= length(trigstruct.allTriggers)
    error('Individual trial types do not sum to number of total trials')
end
