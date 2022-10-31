function fixPosTrackingFiles_cmperpix(topRawDir)

%This function converts data in videoPositionTracking files from pixels to
%cm. This should be done automatically using the range function, but in
%case it isn't this function will fix it and overwrite existing data to cm.
%Requires pixelscale field in videoPositionTracking file with conversion-
%should be added if range function is used. 

cd(topRawDir);
rawDir = dir();
rawDir= {rawDir(3:end).name};

for d = 1:length(rawDir) %loop over days
    dayDir = rawDir{d};
    
    cd(dayDir)

tFiles = dir(['*.videoPositionTracking']);


    for fileInd = 1:length(tFiles) %loop over all position tracking files
        tmpFileName = tFiles(fileInd).name;
        tmpPosData = readTrodesExtractedDataFile(tmpFileName);

        %get pix/cm string
        pixstr = tmpPosData.pixelscale;

        %find the start of the extra characters (not numbers)
        txtstart = strfind(pixstr, 'pix/cm');

        %fix
        pixcm = str2num(str(1:(txtstart-2)));


        for n = 2:length(tmpPosData.fields)
            rawpos = double(tmpPosData.fields(n).data);
            newpos = rawpos/pixcm;
            tmpPosData.fields(n).data = newpos;
        end

    end
end

