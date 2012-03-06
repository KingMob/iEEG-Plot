% ft_group_ieeg_elec_plot(varargin)
%
% Plots one or more subjects' electrodes on a brain (presumably in MNI space)
%
% written by Matthew Davidson, psymatthew@gmail.com

function ft_group_ieeg_elec_plot(varargin)
    % Set defaults
    elec = [];
    elec_files = {};
    stat_files = {};
    data_files = {};
    contrast = [];
    surf_files = {};
    time_windows = {[0 Inf]};
    using_time = false;
    freq_bandwidths = {[0 Inf]};
    using_freqs = false;
    using_analysis_results = false;
    showing_labels = false;
    showing_signif_labels = false;
    new_fig = false;
    fig_title = '';
    output_image = '';
    hemisph = [];
    elec_size_method = '';
    elec_color_method = '';
    elec_edge_color_method = 'signif';
    force_to_side = '';
    force_to_max_x = false;
    viewpoints = {[90 0] [-90 0] [180 -90] [0 90]};
    pval_signif_thr = [];
    pval_max_thr = [];
    foi = [];
    scaling_factor = 1;
    
    parse_inputs(varargin);
    
    
    % Read in electrode info
    if(isempty(elec))
        elec = read_electrode_info(elec_files);
    end
    
    % Alter and remove electrodes
    elec = remove_depth_electrodes(elec);
    if(isempty(force_to_side))
        elec = remove_hemisph_electrodes(elec, hemisph);
    end
    if(isfield(elec, 'pvals') && ~isempty(elec.pvals))
        wh_nan = isnan(elec.pvals);
        elec = remove_electrodes(elec, wh_nan);
    end
    
    if(size(elec.coords, 1) > 0)
        for i=1:length(time_windows)
            for j=1:length(freq_bandwidths)
                if(using_analysis_results)
                    elec = add_analysis_info(elec, stat_files, data_files, using_time, time_windows{i}, using_freqs, freq_bandwidths{j}, foi, pval_signif_thr);
                end
                
                elec = compute_elec_appearances(elec, elec_size_method, elec_color_method, elec_edge_color_method, scaling_factor, pval_signif_thr, pval_max_thr);
                
                if(all(time_windows{1} == [0 Inf]))
                    time_title = '';
                elseif(time_windows{i}(1) ~= time_windows{i}(2))
                    time_title = sprintf(' - %d-%dms', time_windows{i}(:)*1000);
                else
                    time_title = sprintf(' - %dms', time_windows{i}(1)*1000);
                end
                
                if(all(freq_bandwidths{1} == [0 Inf]))
                    freq_title = '';
                elseif(freq_bandwidths{j}(1) ~= freq_bandwidths{j}(2))
                    freq_title = sprintf(' - %d-%d Hz', freq_bandwidths{j}(:));
                else
                    freq_title = sprintf(' - %d Hz', freq_bandwidths{j}(1));
                end
                
                brain_elec_view = plot_brain_and_elecs(surf_files, hemisph, elec, 'title', sprintf('%s%s%s', fig_title, time_title, freq_title), ...
                    'float_to_camera', 5, 'extra_float_to_camera', elec.signif, 'showing_labels', showing_labels, 'showing_signif_labels', showing_signif_labels, ...
                    'force_to_max_x', force_to_max_x, 'force_to_side', force_to_side);
                
                % Save views
                if(~isempty(output_image))
                    [d f e] = fileparts(output_image);
                    mkdir(d);
                    for j=1:length(viewpoints)
                        brain_elec_view(viewpoints{j});
                        fname = fullfile(d, [f sprintf('%s%s - az%d el%d', time_title, freq_title, viewpoints{j}(:)) e]);
                        saveas(gcf, fname);
                    end
                end
            end
        end
    else
        disp('No electrodes found on this hemisphere');
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parse_inputs(args)
        % Parse inputs
        for i=1:length(args)
            if(ischar(args{i}))
                switch(args{i})
                    case 'elec'
                        elec = args{i+1};
                    case 'elec_files'
                        elec_files = cellstr(args{i+1});
                    case 'stat_files'
                        stat_files = cellstr(args{i+1});
                        using_analysis_results = true;
                    case 'data_files'
                        data_files = args{i+1};
                        using_analysis_results = true;
                    case 'contrast'
                        contrast = args{i+1};
                    case {'surf_file' 'surf_files'}
                        surf_files = cellstr(args{i+1});
                    case {'time_windows'}
                        time_windows = args{i+1};
                        using_time = true;
                    case {'freq_bandwidths'}
                        freq_bandwidths = args{i+1};
                        using_freqs = true;
                    case {'showing_labels' 'show_labels' 'labels'}
                        showing_labels = args{i+1};
                    case {'showing_signif_labels' 'show_signif_labels' 'signif_labels'}
                        showing_signif_labels = args{i+1};
                    case 'hemisph'
                        hemisph = args{i+1};
                    case 'new_fig'
                        new_fig = args{i+1};
                    case {'title' 'fig_title'}
                        fig_title = args{i+1};
                    case 'output_image'
                        output_image = args{i+1};
                    case 'elec_size_method'
                        elec_size_method = args{i+1};
                    case 'elec_color_method'
                        elec_color_method = args{i+1};
                    case 'elec_edge_color_method'
                        elec_edge_color_method = args{i+1};
                    case 'viewpoints'
                        viewpoints = args{i+1};
                    case 'pval_signif_thr'
                        pval_signif_thr = args{i+1};
                    case 'pval_max_thr'
                        pval_max_thr = args{i+1};
                    case 'foi'
                        foi = args{i+1};
                    case 'scaling_factor'
                        scaling_factor = args{i+1};
                    case 'force_to_side'
                        force_to_side = args{i+1};
                    case 'force_to_max_x'
                        force_to_max_x = args{i+1};
                end
            end
        end
        
        % Check inputs
        %assert(using_time || using_freqs, 'Have to have some way of sorting electrodes'); 
        assert(all(cellfun(@(x)length(x)==2, time_windows)), 'Time windows param must be a start and a stop in seconds. Inf and -Inf are allowed.');
        assert(all(cellfun(@(x)length(x)==2, freq_bandwidths)), 'Time windows param must be a start and a stop in seconds. Inf and -Inf are allowed.');
        assert(~isempty(surf_files), 'Must supply at least one surf file.');
        assert(length(surf_files) <= 2, 'Cannot supply more than two surf files.');
        assert(~isempty(elec) || ~isempty(elec_files), 'Must supply elec structure or at least one electrode file.');
        if(~isempty(stat_files))
            assert(length(stat_files) == length(elec_files), 'If using, you must supply the same number of stat files as electrode files.');
        end
        if(~isempty(data_files))
            for i=1:length(data_files)
                assert(length(data_files{i}) == length(elec_files), 'If using, you must supply the same number of data files per cell as electrode files.');
            end
        end
        assert(islogical(showing_labels) || isnumeric(showing_labels), 'Show labels param must be true/false.');
        
        % Set hemisphere view
        if(isempty(hemisph))
            if(length(surf_files) == 2)
                hemisph = 'both';
            else
                hemisph = regexp(surf_files, '[r, l]h', 'match');
                hemisph = char(hemisph{:});
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elec = add_analysis_info(elec, stat_files, data_files, using_time, time_window, using_freqs, freq_band, foi, pval_signif_thr)
    persistent prev_data_files prev_datas prev_stat_files prev_stat_data
    
    MIN_SIGNIF_DUR = .02;
    
    using_stat_data = false;
    using_data = false;
    loading_stat_data = false;
    loading_data = false;
    
    if(~isempty(stat_files)), using_stat_data = true; loading_stat_data = true; end
    if(~isempty(data_files)), using_data = true; loading_data = true; end
    determine_data_cache_status();
    determine_stat_cache_status();
    
    num_subjs = max(elec.subj);
    
    for i=1:num_subjs
        if(using_data)
            if(loading_data)
                for j=1:length(data_files)
                    datas{j} = load(data_files{j}{i});
                    names = fieldnames(datas{j});
                    data_field = names{1};
                    datas{j} = datas{j}.(data_field);
                    datas{j}.label = standardize_elec_names(datas{j}.label);
                    
                    cfg = [];
                    cfg.baseline = [-.3 -.05];
                    cfg.baselinetype = 'relchange';
                    datas{j} = ft_freqbaseline(cfg, datas{j});
                end
                prev_datas{i} = datas;
            else
                datas = prev_datas{i};
            end
           
            elec_data_vali = get_data_val(datas, time_window, foi);
        end
        
        if(using_stat_data)
            if(loading_stat_data)
                load(stat_files{i});
                stat.label = standardize_elec_names(stat.label);
                prev_stat_data{i} = stat;
            else
                stat = prev_stat_data{i};
            end
            
            elec_pvalsi = stat.prob;
            elec_stati = stat.stat;
            elec_signifi = stat.mask;
            
            if(using_time)
                % Compute significant duration
                if(time_window(1) ~= time_window(2))
                    wh_window = find(stat.time >= time_window(1) & stat.time <= time_window(2));
                    if(~isempty(pval_signif_thr))
                        signif_dur = sum(stat.prob(:,wh_window) <= pval_signif_thr, 2) / (1/(stat.time(2)-stat.time(1)));
                    else
                        signif_dur = sum(stat.mask(:,wh_window), 2) / (1/(stat.time(2)-stat.time(1)));
                    end
                    elec_signifi = signif_dur > MIN_SIGNIF_DUR;
                else
                    wh_window = nearest(stat.time, time_window(1));
                    if(~isempty(pval_signif_thr))
                        elec_signifi = stat.prob(:,wh_window) <= pval_signif_thr;
                    else
                        elec_signifi = stat.mask(:,wh_window);
                    end
                end
            end
            
            if(using_freqs)
                wh_window = find(stat.freq >= freq_band(1) & stat.freq <= freq_band(2));
                if(~isempty(pval_signif_thr))
                    elec_signifi = any(stat.prob(:,wh_window) <= pval_signif_thr, 2);
                else
                    elec_signifi = any(stat.mask(:,wh_window), 2);
                end
            end
            
            % Compute min p-val
            if(numel(stat.prob) ~= length(stat.label))
                elec_pvalsi = min(stat.prob(:,wh_window), [], 2);
                elec_stati = mean(stat.stat(:,wh_window), 2);
            end
        end
        
        remove_missing_elecs();
        insert_subjs_data();
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function determine_data_cache_status()
        if(~isempty(prev_data_files) && ~isempty(data_files) && length(prev_data_files) == length(data_files))
            for i=1:length(data_files)
                if(length(prev_data_files{i}) ~= length(data_files{i}) || ~all(cellfun(@strcmp, data_files{i}, prev_data_files{i})))
                    return
                end
            end
            loading_data = false;
        end
        prev_data_files = data_files;
    end
    
    function determine_stat_cache_status()
        if(~isempty(prev_stat_files) && ~isempty(stat_files) && length(prev_stat_files) == length(stat_files))
            if(all(cellfun(@strcmp, stat_files, prev_stat_files)))
                loading_stat_data = false;
            end
        end
        prev_stat_files = stat_files;
    end
    
    function remove_missing_elecs()
        if(using_data)
            data_elec_names = datas{1}.label;
            [elec_namesi, wh_valid_elec, wh_valid_data_elec] = intersect(elec.names(elec.subj == i), data_elec_names);
            
            subj_offset = find(elec.subj == i, 1, 'first') - 1;
            wh_invalid_elec = setxor(wh_valid_elec, 1:sum(elec.subj == i));
            elec = remove_electrodes(elec, wh_invalid_elec + subj_offset);

            elec_data_vali = elec_data_vali(wh_valid_data_elec);
        end
        
        if(using_stat_data)
            stat_elec_names = cellstr(char(stat.label));
            [elec_namesi, wh_valid_elec, wh_valid_stat_elec] = intersect(elec.names(elec.subj == i), stat_elec_names); %requires presorting
            
            subj_offset = find(elec.subj == i, 1, 'first') - 1;
            wh_invalid_elec = setxor(wh_valid_elec, 1:sum(elec.subj == i));
            elec = remove_electrodes(elec, wh_invalid_elec + subj_offset);
            
            assert(length(elec.names(elec.subj == i)) == length(stat_elec_names(wh_valid_stat_elec)));
            assert(all(cellfun(@strcmp, elec.names(elec.subj == i), stat_elec_names(wh_valid_stat_elec))));
            
            elec_pvalsi = elec_pvalsi(wh_valid_stat_elec);
            elec_stati = elec_stati(wh_valid_stat_elec);
            elec_signifi = elec_signifi(wh_valid_stat_elec);
        end
    end
    
    function insert_subjs_data()
        if(exist('elec_data_vali', 'var') && ~isempty(elec_data_vali))
            elec.data_val(elec.subj == i) = elec_data_vali;
        end
        if(exist('elec_signifi', 'var') && ~isempty(elec_signifi))
            elec.signif(elec.subj == i) = elec_signifi;
        end
        if(exist('elec_pvalsi', 'var') && ~isempty(elec_pvalsi))
            elec.pvals(elec.subj == i) = elec_pvalsi;
        end
        if(exist('elec_stati', 'var') && ~isempty(elec_stati))
            elec.stat(elec.subj == i) = elec_stati;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function elec_data_val = get_data_val(datas, time_window, foi)
    ft_type = datatype(datas{1});
    wh_window = find(datas{1}.time >= time_window(1) & datas{1}.time <= time_window(2));
    
    switch(ft_type)
        case 'freq'
            wh_freq = (datas{1}.freq == foi);
            switch(length(datas))
                case 2
                    elec_data_val = mean(datas{1}.powspctrm(:, wh_freq, wh_window) - datas{2}.powspctrm(:, wh_freq, wh_window), 2);
                case 1
                    elec_data_val = mean(squeeze(datas{1}.powspctrm(:, wh_freq, wh_window)), 2);
                otherwise
                    error('How many conditions in this thing?');
            end
        case 'timelock'
            switch(length(datas))
                case 2
                    elec_data_val = mean(abs(datas{1}.avg(:, wh_window)) - abs(datas{2}.avg(:, wh_window)), 2);
                case 1
                    disp('Displaying diff from 0');
                otherwise
                    error('How many conditions in this thing?');
            end
        otherwise
            error('Unknown fieldtrip type: %s\n', ft_type);
    end
end

function elec = remove_depth_electrodes(elec)
    wh_depths = wh_str('^D', elec.names);
    elec = remove_electrodes(elec, wh_depths);
end

function elec = remove_hemisph_electrodes(elec, hemisph)
    MEDIAL_DIST = 2;
    
    switch(lower(hemisph))
        case {'l' 'lh'}
            wh_to_remove = elec.coords(:,1) > MEDIAL_DIST;
            elec = remove_electrodes(elec, wh_to_remove);
        case {'r' 'rh'}
            wh_to_remove = elec.coords(:,1) < -MEDIAL_DIST;
            elec = remove_electrodes(elec, wh_to_remove);
        case 'both'
        otherwise
            error('Unknown hemisphere: %s\n', hemisph);
    end
end


function elec = compute_elec_appearances(elec, elec_size_method, elec_color_method, elec_edge_color_method, scaling_factor, pval_signif_thr, pval_max_thr)
    elec = set_elec_size(elec, elec_size_method, scaling_factor);
    elec = set_elec_color(elec, elec_color_method, pval_signif_thr, pval_max_thr);
    elec = set_elec_edge_color(elec, elec_edge_color_method, pval_signif_thr);    
end

function elec = set_elec_size(elec, elec_size_method, scaling_factor)
    DEFAULT_SCALING_FACTOR = 1;
    MAX_RADIUS = 5;
    MIN_RADIUS = 1;
    
    if(~exist('scaling_factor', 'var') || isempty(scaling_factor)), scaling_factor = DEFAULT_SCALING_FACTOR; end
    
    elec.radius = ones(size(elec.names));
    switch(elec_size_method)
        case 'log_signif'
            elec.radius(elec.signif) = -log(elec.pvals(elec.signif));
        case 'log_signif_only_if_signif'
            elec.radius(elec.signif) = -log(elec.pvals(elec.signif));
            elec.radius(~elec.signif) = 0;
        case 'loglog_signif'
            elec.radius(elec.signif) = log(-log(elec.pvals(elec.signif)));
        case 'loglog_signif_only_if_signif'
            elec.radius(elec.signif) = log(-log(elec.pvals(elec.signif)));
            elec.radius(~elec.signif) = 0;
        case 'amplitude'
            elec.radius = abs(elec.data_val);
        case 'log_amplitude'
            elec.radius = log10(abs(elec.data_val));
        case 'amplitude_if_signif'
            elec.radius(elec.signif) = abs(elec.data_val(elec.signif));
        case 'log_amplitude_if_signif'
            elec.radius(elec.signif) = log10(abs(elec.data_val(elec.signif)));
        case 'amplitude_only_if_signif'
            elec.radius(elec.signif) = log10(abs(elec.data_val(elec.signif)));
            elec.radius(~elec.signif) = 0;
        case 'amplitude_if_signif_wmin'
            elec.radius(elec.signif) = max(MIN_RADIUS, log2(abs(elec.data_val(elec.signif)))); % use max to get a min size!
        case 'stat'
            elec.radius = elec.stat / max(elec.stat);
        case 'stat_if_signif'
            elec.radius(elec.signif) = log10(abs(elec.stat(elec.signif)));
    end
    elec.radius = elec.radius * scaling_factor;
    elec.radius(elec.radius > MAX_RADIUS) = MAX_RADIUS;
end

function elec = set_elec_color(elec, elec_color_method, pval_signif_thr, pval_max_thr)
    if(~exist('pval_signif_thr', 'var') || isempty(pval_signif_thr))
        pval_signif_thr = .0001;
    end
    if(~exist('pval_max_thr', 'var') || isempty(pval_max_thr))
        pval_max_thr = .000000001;
    end
        
    elec.color = repmat([0 0 0], length(elec.names), 1);
    switch(elec_color_method)
        case 'data'
            elec.color(elec.data_val > 0, 1) = 1;
            elec.color(elec.data_val < 0, 3) = 1;
        case 'stat'
            elec.color(elec.stat > 0, 1) = 1;
            elec.color(elec.stat < 0, 3) = 1;
        case 'data_if_signif'
            elec.color(elec.signif & elec.data_val > 0, 1) = 1;
            elec.color(elec.signif & elec.data_val < 0, 3) = 1;
        case 'stat_if_signif'
            elec.color(elec.signif & elec.stat > 0, 1) = 1;
            elec.color(elec.signif & elec.stat < 0, 3) = 1;
        case 'signif'
            elec.color = .7 * ones(size(elec.pvals, 1), 3);
            elec.color(elec.signif, :) = repmat([1 1 1], sum(elec.signif), 1);
        case 'pval_tiered'
            hotmap = hot(6);
            elec.color(elec.pvals < .05, :) = repmat(hotmap(1,:), sum(elec.pvals < .05), 1);
            elec.color(elec.pvals < .01, :) = repmat(hotmap(2,:), sum(elec.pvals < .01), 1);
            elec.color(elec.pvals < .001, :) = repmat(hotmap(3,:), sum(elec.pvals < .001), 1);
            elec.color(elec.pvals < .0001, :) = repmat(hotmap(4,:), sum(elec.pvals < .0001), 1);
            elec.color(elec.pvals < .00001, :) = repmat(hotmap(5,:), sum(elec.pvals < .00001), 1);
            elec.color(elec.pvals < .000001, :) = repmat(hotmap(6,:), sum(elec.pvals < .000001), 1);
        case 'signif_pval_tiered'
            p_thrs = [.01 .0001 .000001 .00000001 .0000000001];
            p_colors = {[1 0 0] [1 .6275 0] [1 1 0] [1 1 .5] [1 1 1]};

            for i=1:length(p_thrs)
                elec.color(elec.signif & elec.pvals < p_thrs(i), :) = repmat(p_colors{i}, sum(elec.signif & elec.pvals < p_thrs(i)), 1);
            end
        case 'logpval_saturation'
            MAX_NEGLOG_PVAL = -log(pval_max_thr);
            MIN_NEGLOG_PVAL = -log(pval_signif_thr);
            HUE = 5/6; % magenta
            
            %Set to grey
            elec.color = .7 * ones(size(elec.pvals, 1), 3);
            
            logpvals = -log(elec.pvals);
            wh_to_color = find(logpvals > MIN_NEGLOG_PVAL);
            
            curr_logpval = logpvals(wh_to_color);
            curr_perc = min(max(0, curr_logpval - MIN_NEGLOG_PVAL) / MAX_NEGLOG_PVAL, 1);
            sat = curr_perc;
            val = .7 + curr_perc * .2;
            elec.color(wh_to_color, :) = hsv2rgb([HUE * ones(size(sat)) sat val]);
        case 'logpval_saturation_bidir'
            MAX_NEGLOG_PVAL = -log(pval_max_thr);
            MIN_NEGLOG_PVAL = -log(pval_signif_thr);
            HUES = [0 2/3]; % red, blue
            
            %Set to grey
            elec.color = .7 * ones(size(elec.pvals, 1), 3);
            
            logpvals = -log(elec.pvals);
            wh_to_color = find(logpvals > MIN_NEGLOG_PVAL);

            curr_logpval = logpvals(wh_to_color);
            curr_perc = min(max(0, curr_logpval - MIN_NEGLOG_PVAL) / MAX_NEGLOG_PVAL, 1);
            sat = curr_perc;
            val = .7 + curr_perc * .3;
            hue = HUES((elec.stat(wh_to_color) < 0) + 1)';
            elec.color(wh_to_color, :) = hsv2rgb([hue sat val]);
        case 'logpval_saturation_pos'
            MAX_NEGLOG_PVAL = -log(pval_max_thr);
            MIN_NEGLOG_PVAL = -log(pval_signif_thr);
            HUE = 0; % red
            
            %Set to grey
            elec.color = .7 * ones(size(elec.pvals, 1), 3);
            
            logpvals = -log(elec.pvals);
            wh_to_color = find(logpvals > MIN_NEGLOG_PVAL & elec.stat > 0);
            
            curr_logpval = logpvals(wh_to_color);
            curr_perc = min(max(0, curr_logpval - MIN_NEGLOG_PVAL) / MAX_NEGLOG_PVAL, 1);
            sat = curr_perc;
            val = .7 + curr_perc * .3;
            elec.color(wh_to_color, :) = hsv2rgb([HUE * ones(size(sat)) sat val]);
        case 'logpval_saturation_neg'
            MAX_NEGLOG_PVAL = -log(pval_max_thr);
            MIN_NEGLOG_PVAL = -log(pval_signif_thr);
            HUE = 2/3; % blue
            
            %Set to grey
            elec.color = .7 * ones(size(elec.pvals, 1), 3);
            
            logpvals = -log(elec.pvals);
            wh_to_color = find(logpvals > MIN_NEGLOG_PVAL & elec.stat < 0);
            
            curr_logpval = logpvals(wh_to_color);
            curr_perc = min(max(0, curr_logpval - MIN_NEGLOG_PVAL) / MAX_NEGLOG_PVAL, 1);
            sat = curr_perc;
            val = .7 + curr_perc * .3;
            elec.color(wh_to_color, :) = hsv2rgb([HUE * ones(size(sat)) sat val]);
        case 'none'
        otherwise
            error('Unknown electrode coloring method: %s\n', elec_color_method);
    end
end

function elec = set_elec_edge_color(elec, elec_edge_color_method, pval_signif_thr)
    elec.edge_color = elec.color;
    switch(elec_edge_color_method)
        case 'signif'
            elec.edge_color(elec.signif,:) = 1;
        case 'pval_thr'
            elec.edge_color(elec.pvals < pval_signif_thr,:) = 1;
        case 'none'
    end
end
