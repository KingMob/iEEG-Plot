function newelec = keeponly_electrodes(elec, wh)
    newelec.coords = elec.coords(wh,:);
    newelec.names = elec.names(wh);
    if(isfield(elec, 'pvals')), newelec.pvals = elec.pvals(wh); end
    if(isfield(elec, 'data_val')), newelec.data_val = elec.data_val(wh); end
    if(isfield(elec, 'signif')), newelec.signif = elec.signif(wh); end
    if(isfield(elec, 'stat')), newelec.stat = elec.stat(wh); end
    
    if(isfield(elec, 'subj')), newelec.subj = elec.subj(wh); end
    if(isfield(elec, 'radius')), newelec.radius = elec.radius(wh); end
end