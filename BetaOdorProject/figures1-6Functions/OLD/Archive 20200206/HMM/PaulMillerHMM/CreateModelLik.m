function [P2,Q2]=CreateModelLik(newdat,neurons,taste,trialstouse,time,its,states,diag)
dat = squeeze(newdat(taste,trialstouse,time));

%initialize the Probability matrix uniformly with sojourn values
for i = 1:states
    for j = 1:states
        if i == j
            P(i,j)= diag;
        else
            P(i,j)= (1 - diag)/(states-1);
        end;
    end;
end;

%randomly initialize the firing rate matrix
Q = rand(states, neurons+1);
Q = normalise(Q,1);

save P P; save Q Q; save trialstouse trialstouse;

%call the matlab training function to create the HMM
[P2, Q2, lik] = hmmtrainLog(dat, P, Q,'Verbose', false,'Maxiterations', its,'Algorithm','Baumwelch','Tolerance',1e-6);

save P2 P2;
save Q2 Q2;
save lik lik;
