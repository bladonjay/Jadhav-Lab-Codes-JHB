%This program (quick and dirty) does the following:
%(1) Formats your data to run with my code
%(2) Initializes HMM parameters
%(3) Creates 10 models
%(4) Selects Best Model
%(5) Plots state probabilities for each trial of best model
%LMJ 10/31/2006

%First import your files into matlab
%Name each file as 'n1' for neuron 1, 'n2' etc...
%Once loaded, save all neuron files together >> save neurons;

%*******************
%******************* (1) Format data
%*******************
close all
load neurons; %load all neuron files
neurons = 10; %number of neurons in ensemble
alltastes = 4; %to match other datasets, only really filling in 1 
tottime = 3000; %3 secs
bintime = 15;
%trials = max(max(n1(:,1))); %currently 10 trials
trials = 10;
dat = zeros(neurons,1,trials,tottime/bintime); %preallocate data matrix


%Make a single structure containing file info for each neuron
n(1).data = n1; n(2).data = n2; n(3).data = n3; n(4).data = n4; 
n(5).data = n5; n(6).data = n6; 
n(7).data = n7; n(8).data = n8;
n(9).data = n9; n(10).data = n10;
%n(11).data = n11; n(12).data = n12;

%Parse into trials, digitize spike times into 4-ms bins of 0s,1s
for i = 1:neurons
    for j = 1:trials
        x = ceil(1000/bintime*n(i).data(find(n(i).data(:,1) == j),2));
	if x < tottime+1
           dat(i,1,j,x) = 1;
	end;
    end;
end; 
save dat dat;

%Combine all spikes into 1 array, (some loss of simultaneous spikes)
%where spikes from unit 1 are marked with '1', spikes from unit 2 = '2' ...
for i = 1:neurons
    dat(i,:,:,:) = dat(i,:,:,:)*i;
end; 

if neurons > 1
    for i = 1:alltastes
        for j = 1:trials
            for k = 1:tottime/bintime
                temp = find(dat(:,i,j,k) > 0); 
                if temp
                    temp2 = randperm(length(temp)); 
                    newdat(i,j,k) = squeeze(dat(temp(temp2(1)),i,j,k));
                else
                    newdat(i,j,k) = 0;
                end;
            end;
        end;
    end;
 end; 
 
newdat = newdat + 1; %Data entered into HMMtrain cannot contain zeros 
save newdat newdat;

%*******************
%******************* (2) Initialize HMM parameters
%*******************

%Make HMMs from the data
time = ':'; neurons = ':'; trialstouse = ':'; %use all time, all neurons, all trials
its = 400; %max number of times to run EM algorithm, we use 200-400, play with 100-600
states = 5; %number of states, we use 5, play with 2-10
diag = 0.8; %diagonal of probability matrix (probability of remaining in the same state)
             %we use 0.98, play with 0.90-0.999
%NOTE: the emissions matrix is randomly initialized
taste = 1; %taste will always  = 1;

%*******************
%******************* (3) Create 10 models
%*******************

%We actually make 10 models then choose the most likely model
%to check random initial conditions and that not stuck in local minimum
%Rough estimate: for 8 neurons, 10 trials, 200 its, runtime = 10 mins

mkdir Model_1; cd Model_1; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2; %transition matrix and emissions matrix
load lik; Likes(1) = lik(end); %keep track of model likelihood
cd ..;
 
mkdir Model_2; cd Model_2; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(2) = lik(end);
cd ..;

mkdir Model_3; cd Model_3; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(3) = lik(end);
cd ..;
 
mkdir Model_4; cd Model_4; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(4) = lik(end);
cd ..;
 
mkdir Model_5; cd Model_5; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(5) = lik(end);
cd ..;
 
mkdir Model_6; cd Model_6; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(6) = lik(end);
cd ..;
 
mkdir Model_7; cd Model_7; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(7) = lik(end);
cd ..;
 
mkdir Model_8; cd Model_8; 
[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(8) = lik(end);
cd ..;
 
%mkdir Model_9; cd Model_9; 
%[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
%save P2 P2; save Q2 Q2;
%load lik; Likes(9) = lik(end);
%cd ..;

%mkdir Model_10; cd Model_10; 
%[P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag); 
%save P2 P2; save Q2 Q2;
%load lik; Likes(10) = lik(end);
%cd ..;

%*******************
%******************* (4) Select most likely model
%*******************
[low, modelno] = max(Likes);
bestmod = ['Model_' num2str(modelno)];

%*******************
%******************* (5) Plot Probabilities
%*******************

%For each trial of best model, plot the ensemble rasters 
%overlaid with state probability
load dat; load newdat; 
eval(['cd ', bestmod]);
load P2; load Q2;
eval(['cd ..']);

save bestmod bestmod

[neurons tastes trials time] = size(dat);
time = [1:time]; taste = 1;
for j=1:trials
    figure;
    seq = squeeze(newdat(taste,j,time))';
    PSTATES = hmmdecode(seq, P2, Q2);
    hold on; plot(bintime*time, PSTATES'*neurons);
    savefile = [ 'probs_' num2str(j) '.dat'];
    PSTATES2 = PSTATES';
    eval(['save ' , savefile , ' PSTATES2 -ascii']);
    for k = 1:neurons
        x = squeeze(dat(k,taste,j,time));
        y = find(x == 1); hold on;
        for l = 1:length(y)
            plot(bintime*[y(l) y(l)],[k-1 k-0.5]); 
            hold on;
            axis([-50 3050 -1 neurons + 1]);
        end;
        savefile = [ 'spikes_' num2str(j) '_' num2str(k) '.dat'];
        eval(['save ' , savefile, ' y -ascii']);
    end;
end;        



