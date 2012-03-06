function brain_elec_view = plot_brain_and_elecs(surf_brain, viewpoint, elec, varargin)
    fig_title = '';
    showing_labels = false;
    showing_signif_labels = false;
    orient_to_camera = true;
    force_to_nearest_vertex = true;
    float_to_camera = 5;
    extra_float_to_camera = [];
    force_to_side = '';
    force_to_max_x = false;
    
    parse_inputs(varargin);
    
    if(ischar(elec))
        elec_file = elec;
        elec = [];
        
        [elec.names x y z] = textread(elec_file, '%s %f %f %f');
        elec.coords = [x y z];
    end
    %elec.names{end+1} = 'Origin';
    %elec.coords(end+1,:) = [0 0 0];
    

    viewpoint = check_viewpoint(viewpoint);
    [brain_figh brain_view] = plot_brain(surf_brain, viewpoint);
    elec = alter_electrodes(elec, force_to_nearest_vertex, force_to_side, force_to_max_x); % Needs brain plotted to force to vertices
    update_view_and_elecs(viewpoint);
    
    % Set viewing parameters
    axis tight;
    title(fig_title, 'Color', 'k', 'FontSize', 16);
    brain_elec_view = @update_view_and_elecs;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parse_inputs(args)
        if(~exist('args', 'var') || isempty(args)), return; end
        for i=1:length(args)
            if(ischar(args{i}))
                switch(lower(args{i}))
                    case {'title', 'fig_title'}
                        fig_title = args{i+1};
                    case {'show_labels' 'showing_labels'}
                        showing_labels = args{i+1};
                    case {'show_signif_labels' 'showing_signif_labels'}
                        showing_signif_labels = args{i+1};
                    case 'orient_to_camera'
                        orient_to_camera = args{i+1};
                    case 'force_to_nearest_vertex'
                        force_to_nearest_vertex = args{i+1};
                    case 'force_to_max_x'
                        force_to_max_x = args{i+1};
                    case 'force_to_side'
                        force_to_side = args{i+1};
                    case 'float_to_camera'
                        float_to_camera = args{i+1};
                    case 'extra_float_to_camera'
                        extra_float_to_camera = args{i+1};
                end
            end
        end
    end
    
    function update_view_and_elecs(viewpoint)
        %figure(brain_figh);
        set(0, 'CurrentFigure', brain_figh);
        %delete(findobj(brain_figh, 'Tag', 'electrode'));
        delete(findobj(gca, 'Tag', 'electrode'));
        brain_view(viewpoint);
        plot_electrodes(elec, viewpoint, orient_to_camera, float_to_camera, showing_labels, showing_signif_labels, extra_float_to_camera);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function elec = alter_electrodes(elec, force_to_nearest_vertex, force_to_side, force_to_max_x)
    switch(upper(force_to_side))
        case 'R'
            elec.coords(:,1) = abs(elec.coords(:,1));
        case 'L'
            elec.coords(:,1) = -abs(elec.coords(:,1));
    end
    
    if(force_to_nearest_vertex)
        elec.coords = compute_nearest_vertices(elec.coords);
    end
    
    if(force_to_max_x)
        max_x = max(abs(elec.coords(:,1)));
        elec.coords(:,1) = max_x * sign(elec.coords(:,1));
        elec.coords(elec.signif,1) = elec.coords(elec.signif,1) + sign(elec.coords(elec.signif,1));
    end
end

function viewpoint = check_viewpoint(viewpoint)
    if(ischar(viewpoint))
        switch(upper(viewpoint))
            case {'LH', 'L'}
                viewpoint = [270, 0];
            case {'RH', 'R'}
                viewpoint = [90, 0];
            case 'TOP'
                viewpoint = [180, 90];
            case {'BOTH' 'BOTTOM'}
                viewpoint = [180, -90];
            otherwise
                error('Unknown viewpoint option %s\n', viewpoint);
        end
    end
end


function nearest_vertices = compute_nearest_vertices(coords)
    persistent prev_coords prev_nearest_vertices
    
    if(~exist('prev_coords', 'var') || isempty(prev_coords) || size(coords, 1) ~= size(prev_coords, 1) || ~all(coords(:) == prev_coords(:)))
        vertices = get_fig_vertex_info();
        
        [D wh_closest] = pdist2(vertices, coords, 'euclidean', 'Smallest', 1);
        nearest_vertices = vertices(wh_closest,:);
        
        prev_coords = coords;
        prev_nearest_vertices = nearest_vertices;
    else
        nearest_vertices = prev_nearest_vertices;
    end
end


function [vertices vertex_normals] = get_fig_vertex_info()
    vertices = [];
    vertex_normals = [];
    hbrain = findobj(gca, 'Tag', 'brain');
    for i=1:length(hbrain)
        vertices = [vertices; get(hbrain(i), 'Vertices')];
        vertex_normals = [vertex_normals; get(hbrain(i), 'VertexNormals')];
    end
end