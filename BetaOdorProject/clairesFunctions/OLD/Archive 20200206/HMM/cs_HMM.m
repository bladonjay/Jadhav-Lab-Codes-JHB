%cs_HMM
%from Paul Miller's code

close all
clear
animID = 7;
day = 5;
region = 'CA1';
win = 1; %seconds
bintime = 10;


[dat, hmmdat] = cs_createHMMmatrix(animID, day, region, win, bintime);

time = ':'; numcells = ':'; trialstouse = ':'; %use all time, all neurons, all trials
its = 400; %max number of times to run EM algorithm, we use 200-400, play with 100-600
states = 3; %number of states, we use 5, play with 2-10
diag = 0.98; %diagonal of probability matrix (probability of remaining in the same state)
             %we use 0.98, play with 0.90-0.999
             
num = 1; %num will always  = 1;
hmmdir = fileparts(which('cs_HMM')); cd(hmmdir)
mkdir Model_1; cd Model_1; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2; %transition matrix and emissions matrix
load lik; Likes(1) = lik(end); %keep track of model likelihood
cd ..;
 
mkdir Model_2; cd Model_2; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(2) = lik(end);
cd ..;

mkdir Model_3; cd Model_3; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(3) = lik(end);
cd ..;
 
mkdir Model_4; cd Model_4; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(4) = lik(end);
cd ..;
 
mkdir Model_5; cd Model_5; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(5) = lik(end);
cd ..;
 
mkdir Model_6; cd Model_6; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(6) = lik(end);
cd ..;
 
mkdir Model_7; cd Model_7; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(7) = lik(end);
cd ..;
 
mkdir Model_8; cd Model_8; 
[P2,Q2]=CreateModelLik(hmmdat,numcells,num,trialstouse,time,its,states,diag); 
save P2 P2; save Q2 Q2;
load lik; Likes(8) = lik(end);
cd ..;
 
%mkdir Model_9; cd Model_9; 
%[P2,Q2]=CreateModelLik(hmmdat,neurons,num,trialstouse,time,its,states,diag); 
%save P2 P2; save Q2 Q2;
%load lik; Likes(9) = lik(end);
%cd ..;

%mkdir Model_10; cd Model_10; 
%[P2,Q2]=CreateModelLik(hmmdat,neurons,num,trialstouse,time,its,states,diag); 
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
eval(['cd ', bestmod]);
load P2; load Q2;
eval(['cd ..']);

save bestmod bestmod

[numcells, ~, trials, time] = size(dat);
time = [1:time]; num = 1;
for j=1:trials
    figure;
    seq = squeeze(hmmdat(num,j,time))';
    PSTATES = hmmdecode(seq, P2, Q2);
    hold on; plot(bintime*time, PSTATES'*numcells);
    savefile = [ 'probs_' num2str(j) '.dat'];
    PSTATES2 = PSTATES';
    eval(['save ' , savefile , ' PSTATES2 -ascii']);
    for k = 1:numcells
        x = squeeze(dat(k,num,j,time));
        y = find(x ~= 0); hold on;
        for l = 1:length(y)
            plot(bintime*[y(l) y(l)],[k-1 k-0.5]); 
            hold on;
            axis([-50 (win*1000)+50 -1 numcells + 1]);
        end
        savefile = [ 'pstates_' num2str(j) '_' num2str(k)];
        save(savefile,'PSTATES');
        %eval(['save ' , savefile, ' y -ascii']);
    end
end       