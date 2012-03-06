function elec = reorder_electrodes(elec, subj, wh)
    if(length(wh) ~= sum(elec.subj == subj))
        error('Unable to reorder electrodes: wh var is wrong size');
    end
    
    wh_subj = elec.subj == subj;
    subj_offset = find(elec.subj == subj, 1, 'first') - 1;
    wh = wh + subj_offset;
    
    elec.coords(wh_subj,:) = elec.coords(wh,:);
    elec.names(wh_subj) = elec.names(wh);
    if(isfield(elec, 'pvals')), elec.pvals(wh_subj) = elec.pvals(wh); end
    if(isfield(elec, 'data_val')), elec.data_val(wh_subj) = elec.data_val(wh); end
    if(isfield(elec, 'signif')), elec.signif(wh_subj) = elec.signif(wh); end
    if(isfield(elec, 'stat')), elec.stat(wh_subj) = elec.stat(wh); end
    
    if(isfield(elec, 'radius')), elec.radius(wh_subj) = elec.radius(wh); end
end
