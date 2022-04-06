function epochs = cs_getRunEpochs(animDir, animal, runtype)
%runtype = odorplace;
taskfiles = dir([animDir,animal,'task*']);

epochs = [];
for f = 1:length(taskfiles)
    load([animDir,taskfiles(f).name]);
    
    runfilter = ['strcmp($environment,''',runtype,''')'];
    ep = evaluatefilter(task, runfilter);
    epochs = [epochs;ep];
end
end
