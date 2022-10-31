function [params] = sjcs_baselinespecgram(prefix, days, epochs, tets, do_wrtgnd,varargin)
% Shantanu - Nov 2012. 
% From Kenny and Maggies - event_spectrograms .m and calcriptriggerredspectrograms.m respectively
% This is for getting baseline values for given eeg tets for normalization when needed. Save These. 

%CS - cleaned it up, edited for my data 11/17/17

if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=1; %% Epochs 
end
if nargin<4,
    tets=1; %
end
if nargin<5,
    do_wrtgnd=0; % Whether to also do with respect to ground
end


% Define Chronux params
% -------------------------------------------
%movingwin = [100 10]/1000;                
params.Fs = 1500;
% params.fpass = [0 40]; % params.fpass = [0 400];
params.tapers = [1 1];
params.err = [2 0.05];

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'movingwin'
            movingwin = varargin{option+1};
        case 'fpass'
    	    params.fpass = varargin{option+1};
    end
end

savetag = 'tmp';

if params.fpass(2) == 400
    savetag = '';
    movingwin = [100 10]/1000; 
end
if params.fpass(2) == 100
    savetag = 'mid';
    movingwin = [400 20]/1000;
end
if params.fpass(2) == 40
    movingwin = [1000 20]/1000;
    savetag = 'low';
end
if params.fpass(2) <= 20
    movingwin = [1000 20]/1000;
    savetag = 'floor';
end

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'movingwin'
            movingwin = varargin{option+1};
        case 'fpass'
    	    params.fpass = varargin{option+1};
    end
end

% SET DATA
% -------------------------------------------

[topDir]= cs_setPaths();
directoryname = [topDir, prefix, 'Expt\', prefix, '_direct\'];


% Get EEGs and spectrogram it
% ------------------------------
cd([directoryname,'\EEG\']);
savedir = [directoryname,'\EEGSpec\'];

for d=1:length(days)
    
    day = days(d);
    daystring = getTwoDigitNumber(day);
    
    for t=1:length(tets)
        tet=tets(t);
        
        tetstring = getTwoDigitNumber(tet);
        
        eegspec = []; dummy=[]; % Save 1 file for each tet in a day 
        if do_wrtgnd==1
            eeggndspec = []; dummy_gnd = [];
        end
         
        flaggnd=0;
        for e=1:length(epochs)
            ep = epochs(e);
            
            epstring = getTwoDigitNumber(ep);
            
            eeg=[]; eeggnd=[];
            epoch=epochs(e);
%                      
            if do_wrtgnd==1
                curreeggndfile = [directoryname,'\EEG\',prefix,'eeg',daystring,'-',epstring,'-',tetstring];
                if (exist([curreeggndfile,'.mat'],'file'))==2
                    disp(['Doing Day ',num2str(day) ', Ep ',num2str(epoch),', Tet ',num2str(tet)]);
                    load(curreeggndfile);
                    specgram_gnd = mtspecgramc(eeg{day}{epoch}{tet}.data,movingwin,params);
%                     eeggndspec{day}{epoch}{tet}.meanspec = mean(specgram_gnd,1);
%                     eeggndspec{day}{epoch}{tet}.stdspec = std(specgram_gnd,1);
                    dummy_gnd=[dummy_gnd;specgram_gnd];
                    % Dont save specgram in file - its too big
                    %eeggndspec{day}{epoch}{tet}.specgram  = [];
                end
            else
                if exist([prefix,'eegref',daystring,'-',epstring,'-',tetstring,'.mat']) == 2
                disp(['Doing Day ',num2str(day) ', Ep ',num2str(epoch),', Tet ',num2str(tet)]);
                %curreegfile = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',tetstring];
                curreegfile = [directoryname,'\EEG\',prefix,'eegref',daystring,'-',epstring,'-',tetstring,'.mat'];
                load(curreegfile);
                %eegspec{day}{epoch}{tet}.specgram = mtspecgramc(eeg{day}{epoch}{tet}.data,movingwin,params);
                specgram = mtspecgramc(eegref{day}{epoch}{tet}.data,movingwin,params);
                % Dont save specgram in file - its too big
                % --------------------------------------------
%                 eegrefspec{day}{epoch}{tet}.meanspec = mean(specgram,1);
%                 eegrefspec{day}{epoch}{tet}.stdspec = std(specgram,1);
                dummy=[dummy; specgram]; % For combining across epochs 
                end
            % Dont save specgram in file - its too big
            % --------------------------------------------
            %eegspec{day}{epoch}{tet}.specgram  = [];
            end
        end % end epochs
        % Also save mean and std for whole day. Save in fields for 1st epoch - [FIXED This has become last epoch now - by accident. Fix later]
        if do_wrtgnd==1  
            
            eeggndspec{day}{1}{tet}.meandayspec=mean(dummy_gnd,1);
            eeggndspec{day}{1}{tet}.stddayspec=std(dummy_gnd,1);
            % Save File for current day and tet
            savefile = [savedir,prefix,'eeggndspec',savetag,daystring,'-',tetstring];
            save(savefile,'eeggndspec');
            clear eeggndspec
        else
            
            eegrefspec{day}{1}{tet}.meandayspec=mean(dummy,1);
            eegrefspec{day}{1}{tet}.stddayspec=std(dummy,1);
            % Save File for current day and tet
            savefile = [savedir,prefix,'eegrefspec',savetag,daystring,'-',tetstring];
            save(savefile,'eegrefspec');
            clear eegrefspec
        end
        
        
    end % end tets
    

    
end % end days


