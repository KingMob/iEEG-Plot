function elec = read_electrode_info(elec_files)
    if(ischar(elec_files))
        elec_files = cellstr(elec_files);
    end
    
    elec.names = [];
    elec.coords = [];
    elec.subj = [];
    elec.pvals = [];
    elec.stat = [];
    elec.data_val = [];
    elec.signif = logical([]);
    
    for i=1:length(elec_files)
        [elec_namesi x y z] = textread(elec_files{i}, '%s %f %f %f %*[^\n\r]');
        elec_namesi = standardize_elec_names(elec_namesi);
        
        elec.names = [elec.names; elec_namesi];
        elec.coords = [elec.coords; [x y z]];
        elec.subj = [elec.subj; i * ones(size(elec_namesi))];
        elec.data_val = [elec.data_val; zeros(size(elec_namesi))];
        elec.signif = [elec.signif; false(size(elec_namesi))];
        elec.pvals = [elec.pvals; ones(size(elec_namesi))];
        elec.stat = [elec.stat; zeros(size(elec_namesi))];
        
        [~, sort_idx] = sort(elec_namesi);
        elec = reorder_electrodes(elec, i, sort_idx);
    end
end