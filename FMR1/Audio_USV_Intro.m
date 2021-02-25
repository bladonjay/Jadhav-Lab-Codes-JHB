%% wav file analysis
[fileToRead1,dir1]=uigetfile;
newData1 = importdata(fullfile(dir1,fileToRead1));

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
%%

timewin=[60 120];

window=timewin*fs;


data2=(data(window(1):window(2)));
T = 1/fs;             % Sampling period       
L = length(data2);             % Length of signal
t = (0:L-1)*T;   

y=fft(data2(:,1));

p2=abs(y/L); %
p1=p2(1:L/2+1); % nyquist?

f=fs*(0:(L/2))/L;
figure;

plot(f,p1);
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% using gpu code to get a cwt on some data


[s,f,t]=spectrogram(data,64,[],[35000:200:60000],fs);

figure; imagesc(t,f,abs(s))


figure; imagesc(t,f,SmoothMat2(zscore(abs(s),[],2)>5,[5 5],2));


% zooming in to a putative call

window = [51 51.5];

wininds=round(window.*fs);

[s,f,t]=spectrogram(data(wininds(1):wininds(2)),36,24,[35000:100:60000],fs);
figure; imagesc(t,f,SmoothMat2(abs(s),[10 10],1));
