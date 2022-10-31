function epochs = cs_getRunEpochs(animDir, animal, runtype,day)
%runtype = odorplace;


if (nargin < 4)
    taskfiles = dir([animDir,animal,'task*']);
else 
    daystr = getTwoDigitNumber(day);
    taskfiles = dir([animDir,animal,'task',daystr,'.mat']);
end

epochs = [];
for f = 1:length(taskfiles)
    load([animDir,taskfiles(f).name]);
    
    runfilter = ['strcmp($environment,''',runtype,''')'];
    ep = evaluatefilter(task, runfilter);
    epochs = [epochs;ep];
end
end
