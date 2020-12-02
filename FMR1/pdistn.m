function [p05, p95, pmid, pmode, pmatrix] = pdistn(q, s, muone, background_prob)
% function pdistn calculates the setwise distance between observations

pmatrix = [];

for ov = 1:size(q,2)
    
    qq = q(ov);
    ss = s(ov);
    
    dels=1e-4;
    
    pr  = dels:dels:1-dels;
    
    fac   = log( pr./(1-pr)/exp(muone) ) - qq;
    fac   = exp(-fac.^2/2/ss);
    pd    = dels*(sqrt(1/2/pi/ss) * 1./(pr.*(1-pr)).* fac);
    
    % sumpd = cumsum(pd);
    sumpd = cumtrapz(pd);
    
    lowlimit  = find(sumpd>0.05);
    if(~isempty(lowlimit) )
        lowlimit  = lowlimit(1);
    else
        lowlimit  = 1;
    end
    
    highlimit = find(sumpd>0.95);
    % highlimit = find(sumpd>0.995);
    if(~isempty(highlimit) )
        if(length(highlimit)>1)
            highlimit = highlimit(1)-1;
        else
            highlimit =  highlimit(1);
        end
    else
        highlimit = length(pr);
    end
    
    middlimit = find(sumpd>0.5);
    if(~isempty(middlimit))
        middlimit = middlimit(1);
    else
        middlimit = length(pr);
    end
    
    
    p05(ov)   = pr(lowlimit(1));
    p95(ov)   = pr(highlimit(1));
    pmid(ov)  = pr(middlimit(1));
    [y,i]     = max(pd);
    pmode(ov) = pr(i);
    
    pmatrix =[pmatrix; sumpd];
    
    
end

inte = fix(background_prob/dels);

pmatrix = pmatrix(:, inte);
end
