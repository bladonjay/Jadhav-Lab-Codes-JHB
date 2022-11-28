function recreateTrodesComments()

% used to recreate TrodesComments files for SJ5
% initiate in correct day folder
% creates .trodesComments if one is not present for a specific .rec entry
% reads timestamps from each .videoTimeStamps file

% automatically specifies to use offsets for 2:N .rec entries
% by inserting time reset command


allRecFiles=dir('*.h264');
[~,recInd]= sort({allRecFiles.date});
allRecFiles= allRecFiles(recInd);

for i=1:length({allRecFiles.name})
   recName=strtok(allRecFiles(i).name,'.');
   if isempty(dir([recName '.trodesComments']))
     
     fid = fopen( [recName '.trodesComments'],'w');
  
     timeStamps= readCameraModuleTimeStamps([recName '.1.videoTimeStamps']);
     startEpoch= timeStamps(1)*30000;
     endEpoch= timeStamps(end)*30000;
     
     %if i>1;
     %fprintf(fid, '%s\n','time reset')
     %end
     
     fprintf(fid, '%u %s\n',startEpoch,'epoch start')
     fprintf(fid, '%u %s\n',endEpoch,'epoch end')


     fclose(fid);
     
   end
end
    
end
