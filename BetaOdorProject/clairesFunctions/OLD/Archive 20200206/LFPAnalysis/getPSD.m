function [PSD,f,t,E] = getPSD(X,fs,varargin)
    % [PSD,f] = getPSD(X,fs,NAME,VALUE)
    % [PSD,f,t] = getPSD(____,'specgram',1,'psdType','cwt') % time is window centers
    % breaks a signal into windows, computes the power spectral density for
    % each window and averages them returning the resulting spectrum. By
    % default this uses the fourier transform ('fft') with 2 sec windows with
    % 50% overlap with each window being multiplied by the hanning function. 
    % NAME-VALUE pairs
    %   'psdType'   - determines type of transform to get PSD, can be
    %                 'fft','welch','stft','chronux','hilbert' or 'cwt' (default 'stft')
    %
    %   'winSize'   - is the size of each window in sec (default 2)
    %
    %   'noverlap'  - is the amount of overlap in sec (default 1)
    %
    %   'windowFun' - is a function handle of the windowing function to apply to each window (default @hann)
    %                 this is only applied with fft
    %
    %   'freqRange' - array defining the freqeuncy range to focus on (default [0.5 55])
    %
    %   'freqRes'   - the size of frequency steps to return (default 0.25): final
    %                 frequency is interpolated from the output of the PSD
    %                 function. If set to [] then interpolation is not done &
    %                 freqRange won't be enforced for fft.
    %
    %   'specgram'  - flag for whether to retrurn a spectrogram or not 
    %                 (default 0). Only works with fft, cwt, stft or chronux.
    %
    %   'bootstrap' - 1 if you want the error on the spectrum to be
    %   obtained via bootstrap. (default=0, no error computed).
    %   'nboot'     - number of times to bootstrap. (default=10000)
    
    %initialize variables
    psdType = 'stft';
    winSize = 2;
    noverlap = 1;
    windowFun = @hann;
    freqRange = [0.5 55];
    freqRes = 0.25;
    specgram = 0;
    nboot = 10000;
    bootstrapFlg = 0;

    for option = 1:2:length(varargin)-1       
        switch(varargin{option})
            case 'psdType' 
                psdType = varargin{option+1};
        end
    end
    %assignVars(varargin{:});

    N  = numel(X);
    t = winSize/2:winSize-noverlap:(N/fs-winSize/2);
    winSize = fix(winSize*fs);
    noverlap = fix(noverlap*fs);

    % Add helper variables
    nStep = winSize-noverlap;
    f = freqRange(1):freqRes:freqRange(2);

    switch lower(psdType)
        case 'fft'
            [winX,winIdx] = slidingWindow(X,winSize,nStep);
            t = (sum(winIdx)-1)./fs;
            W = windowFun(winSize);
            winX = winX.*W;
            Px = fft(winX,N);
            mx = abs(Px).^2;
            mx = mx./(W'*W);
            npts = N/2+1;
            mx = mx(1:npts,:);
            mx(2:end-1,:) = mx(2:end-1,:).*2;
            Pxx = mx./fs;
            Pxx = Pxx./2;
            if ~specgram
                Pxx = mean(Pxx,2);
            end
            Fx = (0:npts-1)*fs/N;
        case 'welch'
            Pxx = pwelch(X,winSize,noverlap,f,fs);
            Fx = f;
        case {'stft','chronux'}
            params = struct('tapers',[3 4],'Fs',fs,'fpass',freqRange);
            [Pxx,t,Fx] = mtspecgramc(X,[winSize/fs nStep/fs],params);
            if ~specgram && size(Pxx,1)>1
                if bootstrapFlg
                    bootDat = bootstrapDat(Pxx,nboot);
                    Pxx = [bootDat.mean];
                    Pxx_Err = [bootDat.SEM];
                else
                    Pxx = mean(Pxx);
                end
            elseif size(Pxx,1)>1
                Pxx = Pxx';
            end
        case 'cwt'
            [WT,Fx] = cwt(X,fs,'VoicesPerOctave',4);
            Fx = Fx(end:-1:1);
            WT = WT(end:-1:1,:);
            Pxx = abs(WT.*conj(WT));
            clear WT
            t = 0:1/fs:(N-1)/fs;
            if ~specgram
                Pxx = mean(Pxx,2);
            end
        case 'hilbert'
            fBands = [(f-freqRes/2)' (f+freqRes/2)'];
            Pxx = zeros(numel(f),numel(t));
            for k=1:numel(f)
                [b,a] = butter(3,fBands(k,:)/(fs/2),'bandpass');
                fX = filtfilt(b,a,X);
                if all(isnan(fX)) || all(fX==0)
                    error('filter for hilbert transform failed')
                end
                eX = envelope(fX);
                winX = slidingWindow(eX,winSize,nStep);
                Pxx(k,:) = abs(mean(winX).^2);
            end
            if ~specgram
                Pxx = mean(Pxx,2);
            end
            Fx = f;
        otherwise
            PSD = [];
            f = [];
            disp('Sorry that psdType is not available.')
            return;
    end
    PSD = interp1(Fx,Pxx,f,'linear','extrap');
    if ~specgram && exist('Pxx_Err','var')
        Pxx_Err = interp1(Fx,Pxx_Err,f,'linear','extrap');
        t = Pxx_Err;
        
    end
