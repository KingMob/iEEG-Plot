function elec = remove_electrodes(elec, wh)
    elec.coords(wh,:) = [];
    elec.names(wh) = [];
    if(isfield(elec, 'pvals')), elec.pvals(wh) = []; end
    if(isfield(elec, 'data_val')), elec.data_val(wh) = []; end
    if(isfield(elec, 'signif')), elec.signif(wh) = []; end
    if(isfield(elec, 'stat')), elec.stat(wh) = []; end
    
    if(isfield(elec, 'subj')), elec.subj(wh) = []; end
    if(isfield(elec, 'radius')), elec.radius(wh) = []; end
end
