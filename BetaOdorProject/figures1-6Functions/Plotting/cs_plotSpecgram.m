function [maxval, minval] = cs_plotSpecgram(eventTrigSpecgram_file, varargin)

%cd('E:\AnalysesAcrossAnimals\'); %home computer
%cd('D:\OdorPlaceAssociation\AnalysesAcrossAnimals\'); %lab computer

load(eventTrigSpecgram_file);

%----- Params ----- %
region = eventTrigSpecgramData.region;
trigtypes = eventTrigSpecgramData.trigtypes;
freqband = eventTrigSpecgramData.freqband;
    switch freqband
        case 'low'
            freqs = [1:40];
        case 'mid'
            freqs = [1:100];
        case 'floor'
            freqs = [1:15];
    end
    

win = eventTrigSpecgramData.win;
tapers = eventTrigSpecgramData.tapers;
movingwin = eventTrigSpecgramData.movingwin;
reforgnd = eventTrigSpecgramData.ReforGnd;



timewin = (win(1) - (-win(2)));

strstart = (strfind(eventTrigSpecgram_file, 'eventTrigSpecgramData') + length('eventTrigSpecgramData_'));
paramstring = eventTrigSpecgram_file(strstart:(strfind(eventTrigSpecgram_file, '.mat')-1));

trigstring = '';

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'Triggers'
            Triggers = varargin{option+1};
            if strcmp(Triggers,'odorTriggers')
                 trigstring = '';
             else
                 trigstring = [Triggers,'_'];
            end
    end
end


%----- Calculate Averages -----%

std=4;
for g = 1:length(trigtypes) %for each trigtype
    
    
    trigtype = trigtypes{g};
    
    specdata = eventTrigSpecgramData.(trigtype)';
    times = [-win(1): timewin/size(specdata,2): win(2)-(timewin/length(specdata))];
%     binsize = times(2)-times(1);
%     bandwidth = freqs(end)/size(specdata,1);
%     freqs = interp1(freqs, freqs, 1:bandwidth:freqs(end));
%     if length(freqs) < size(specdata,1)
%         freqs(end+1) = freqs(end) + bandwidth;
%     end
    
    %smoothing
    %smoothedspecdata = interp2(times,freqs, specdata, times(1):binsize/2:times(end), freqs(1):bandwidth/2:freqs(end));
    s = gaussian2(std,(2*std));
    smoothedspecdata = filter2(s,specdata); 
    %smoothedspecdata = specdata;
    
    newlength = size(smoothedspecdata,1) - 2*tapers(2);  % size to keep
    newwidth = size(smoothedspecdata,2) - 2*tapers(2);
    smoothedspecdata = smoothedspecdata(tapers(2):newlength, tapers(2):newwidth);
%     times = times(tapers(2):newwidth);
%     freqs = freqs(tapers(2):newwidth);
    
    figure, hold on
    colormap(jet);
    imagesc(times, freqs,smoothedspecdata)
    set(gca,'YDir','normal')
    plot([0 0],[1 freqs(end)],'k--', 'LineWidth', 1.5);
    axis([-win(1) win(2) 1 freqs(end)])
    colorbar
    
%     figure, colormap(jet);
%     imagesc(times, freqs, specdata)
%     set(gca,'YDir','normal')
%     colorbar
%     hold on
%     plot([0 0],[1 40],'k--', 'LineWidth', 1.5);
%     axis([-win(1) win(2) 1 40])

    figtitle1 = ['Spectrogram_',paramstring];
    title(sprintf('%s%d',figtitle1), 'Interpreter', 'none'); %allows underscores in fig titles
    
    [~,figdir] = cs_setPaths();
    figdir = [figdir,'Specgrams\'];
    
    figfile = [figdir,figtitle1];
    
    %saveas(gcf,figfile,'fig');
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    
    maxval = max(smoothedspecdata(:));
    minval = min(smoothedspecdata(:));
    
    clear specdata
end