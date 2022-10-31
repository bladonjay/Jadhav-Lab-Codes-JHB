function out = cs_calcPhaseLocking(sph)

if length(sph) > 10
    stats = rayleigh_test(sph); 
    [moddepth, peakphase] = modulation(sph);
    peakphase_deg = peakphase*(180/pi);
    
    % Von Mises Distribution - From Circular Stats toolbox
    [prefdir, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
    prefdir_deg = prefdir*(180/pi);
    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
else
    sph = [];
    stats=0;
    moddepth=0;
    peakphase =NaN;
    peakphase_deg = NaN;
    kappa=0;
    prefdir = NaN;
    prefdir_deg=NaN;
    prayl=NaN;
    zrayl=0; 
    
end

out.sph = sph;
out.stats = stats;
out.moddepth = moddepth;
out.peakphase = peakphase;
out.peakphase_deg = peakphase_deg;
out.kappa = kappa;
out.prefdir = prefdir;
out.prefdir_deg = prefdir_deg;
out.prayl = prayl;
out.zrayl = zrayl;

end