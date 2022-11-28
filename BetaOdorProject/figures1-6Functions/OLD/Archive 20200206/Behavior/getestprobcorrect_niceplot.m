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

function [pc, lt] = getestprobcorrect_niceplot(bp, background_prob, startflag,plotfig)


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

    [p, q, s, qold, sold] = recfilter(I, sige, qguess, sigsqguess, muone);

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
        if( a1 < cvgce_crit & a2 < cvgce_crit & startflag >= 1)
            fprintf(2, 'EM estimates of RW variance and start point converged after %d steps   \n',  jk)
            break
        elseif ( a1 < cvgce_crit & startflag == 0)
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

if plotfig ==1 
    
t=1:size(p,2)-1;

figure;  clf;
set(gcf,'Position',[100 100 1250 895]);
%CS31-NovelOdors
%cumtr=[1,16,42,76,109,137,166,198,227,264,304,341,373];

% CS15 - RuleSwitch
%cumtr=[1,59,119,167,208,247,278,310,344,378,416];

%CS33 - PreTraining
cumtr = [1,36, 88, 155,226,284,370,411];

oddsess = (1:2:length(cumtr));
evensess = (2:2:length(cumtr));
% for i = 1:(length(cumtr)-1)
%     if any(i == oddsess)
%     patch([cumtr(i), cumtr(i), (cumtr(i)+cumtr(i+1))/2, cumtr(i+1), cumtr(i+1)], ...
%         [0, 1, 1.05, 1, 0], [204,255,229]/256, 'EdgeColor','None'); 
%     else
%     patch([cumtr(i), cumtr(i), (cumtr(i)+cumtr(i+1))/2, cumtr(i+1), cumtr(i+1)], ...
%         [0, 1, 1.05, 1, 0], [204,204,255]/256, 'EdgeColor','None'); 
%     end
% end
    
% patch([1, 1, (16+1)/2, 16, 16], [0, 1, 1.05, 1, 0], [204,255,229]/256, 'EdgeColor','None'); 
% patch([16, 16, (16+40)/2, 40, 40], [0, 1, 1.05, 1, 0], [204,204,255]/256, 'EdgeColor','None'); 
% patch([40, 40, (40+72)/2, 72, 72], [0, 1, 1.05, 1, 0], [204,255,229]/256, 'EdgeColor','None'); 
% patch([72, 72, (72+100)/2, 100, 100], [0, 1, 1.05, 1, 0], [204,204,255]/256, 'EdgeColor','None'); 
% patch([100, 100, (100+132)/2, 132, 132], [0, 1, 1.05, 1, 0], [204,255,229]/256, 'EdgeColor','None'); 
% patch([132, 132, (132+161)/2, 161, 161], [0, 1, 1.05, 1, 0], [204,204,255]/256, 'EdgeColor','None'); 
% patch([161, 161, (161+198)/2, 198, 198], [0, 1, 1.05, 1, 0], [204,255,229]/256, 'EdgeColor','None'); 
% patch([326, 326, (198+361)/2, 361, 361], [0, 1, 1.05, 1, 0], [204,204,255]/256, 'EdgeColor','None'); 


patch([t fliplr(t)], [b05(2:end),fliplr(b95(2:end))],[.7 .7 .7], 'EdgeColor','None')
hold on;

plot(t, bmode(2:end),'k-', 'LineWidth',3);

plot(t, b05(2:end),'Color',[.7 .7 .7]);
plot(t, b95(2:end), 'Color',[.7 .7 .7]);
hold on; [y, x] = find(bp > 0);

%fillspace = [b05(2:end); b95(2:end) ];
%tt = [t; t];
%patch(tt, fillspace, 'b');

box off
% h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
% set(h, 'MarkerEdgeColor', 'k');
% hold on; [y, x] = find(bp == 0);
% h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
% set(h, 'MarkerEdgeColor', 'k');
axis([1 t(end)  0 1.05]);
plot([1 t(end)], [background_prob  background_prob ],'k--')
%line([1 t(end)], [background_prob  background_prob ]);
%title(['IO(0.95) Learning Trial = ' num2str(lt) ' RW variance = ' num2str(sige^2) ]);
xlabel('Trial Number','FontSize',28)
ylabel('Probability of a Correct Response', 'FontSize',28)
title('');
set(0,'defaultaxesfontsize',28);

% xvals = [x,x];        % repeat x values
% yvals = [yy1,yy2];   % vector of upper & lower boundaries
% fill(x,yy,'b') 


% 
% subplot(212)
% plot(t,1 - pmatrix(2:end),'k')
% line([ 1 t(end)],[0.90 0.90]);
% line([ 1 t(end)],[0.99 0.99]);
% line([ 1 t(end)],[0.95 0.95]);
% axis([1 t(end)  0 1]);
% grid on;
% xlabel('Trial Number')
% ylabel('Certainty')

end
pc = [bmode' b05' b95'];
