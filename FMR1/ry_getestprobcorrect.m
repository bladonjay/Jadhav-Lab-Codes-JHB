% [probcorrect] = getestprobcorrect(behavperform, backgroundprob, initialcond)
%       goes through the ones and zeros of behavperform and returns the
%       estimate, with confidence bounds, of the probability of a correct
%	response at each trial using the Smith et al. (2003) algorithm with
%	the background probability of correct.
%
%	initialcond = 0 specifies that the initial condition is fixed at the
%			background probability
%		      1 specifies that the initial condition should be
%			estimated
%		      2 removes xo from likelihood; anposcorrow strong bias
%   varargin: backest = 0 if only want estimate to move forward
%                     = 1 if want to estimate in forward and backward direction
%                     default is 1
%
%   probcorrect
%       The first column of probcorrect is the mode
%       The second column of probcorrect is the lower 5% confidence bound
%       The third column of probcorrect is the upper 5% confidence bound
%       NOTE: first row should be exluded

function [pc, lt] = ry_getestprobcorrect(bp, background_prob, startflag)


% set the total number that could be correct at each trial
nposcorr = ones(size(bp));

I          = [bp' ; nposcorr'];

%starting guess for sige = sqrt(sigma_eps squared)
sige        = sqrt(0.005) ;
sigsqguess  = sige^2;

%set the value of mu from the chance of correct
muone           = log(background_prob/(1-background_prob)) ;

%convergence criterion for sigma_eps_squared
cvgce_crit = 1e-8;

%----------------------------------------------------------------------------------
%loop through EM algorithm

qguess         = 0;  %starting point for random walk q
number_steps  =  2000;

for jk=1:number_steps

    %do forward filter then backward filter, then EM

    [p, q, s, qold, sold] = ry_recfilter(I, sige, qguess, sigsqguess, muone);

    [qnew, signewsq, a]   = backest(q, qold, s, sold);
    
    if (startflag == 1)
        qnew(1) = 0.5*qnew(2);   %updates the initial value of the latent process
        signewsq(1) = sige^2;
    elseif(startflag == 0)
        qnew(1) = 0;             %fixes initial value (no bias at anposcorr)
        signewsq(1) = sige^2;
    elseif(startflag == 2)
        qnew(1) = qnew(2);       %xo = x1 means no prior chance probability
        signewsq(1) = signewsq(2);
    end


    [newsigsq(jk)]         = em_bino(I, qnew, signewsq, a, muone, startflag);


    qnew1save(jk) = qnew(1);

    %check for convergence
    if(jk>1)
        a1 = abs(newsigsq(jk) - newsigsq(jk-1));
        a2 = abs(qnew1save(jk) -qnew1save(jk-1));
        if( a1 < cvgce_crit && a2 < cvgce_crit && startflag >= 1)
            fprintf(2, 'EM estimates of RW variance and start point converged after %d steps   \n',  jk)
            break
        elseif ( a1 < cvgce_crit && startflag == 0)
            fprintf(2, 'EM estimate of RW variance converged after %d steps   \n',  jk)
            break
        end
    end

    sige   = sqrt(newsigsq(jk));
    qguess = qnew(1);
    sigsqguess = signewsq(1);

end


if(jk == number_steps)
    fprintf(2,'failed to converge after %d steps; convergence criterion was %f \n', jk, cvgce_crit)
end

%-----------------------------------------------------------------------------------
%integrate and do change of variables to get confidence limits

[b05, b95, bmid, bmode, pmatrix] = pdistn(qnew, signewsq, muone, background_prob);

%-------------------------------------------------------------------------------------
%find the last point where the 90 interval crosses chance
%for the backward filter (lt)

lt = find(b05 < background_prob);

if(~isempty(lt))
    if(lt(end) < size(I,2) )
        lt = lt(end);
    else
        lt = NaN;
    end
else
    lt = NaN;
end

%-------------------------------------------------------------------------------------
%plot up the figures
% 
t=1:size(p,2)-1;


pc = [bmode' b05' b95'];
