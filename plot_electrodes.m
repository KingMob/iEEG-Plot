function elec_view = plot_electrodes(elec, viewpoint, orient_to_camera, float_to_camera, showing_labels, showing_signif_labels, extra_float_to_camera)
    DEFAULT_RADIUS = 1.5;
    DEFAULT_COLOR = [1 1 1];
    DEFAULT_EDGE_COLOR = [0 0 0];
    
    if(~isfield(elec, 'radius') || isempty(elec.radius))
        elec.radius = DEFAULT_RADIUS * ones(size(elec.names));
    end
    if(~isfield(elec, 'color') || isempty(elec.color))
        elec.color = repmat(DEFAULT_COLOR, length(elec.names), 1);
    end
    if(~isfield(elec, 'edge_color') || isempty(elec.edge_color))
        elec.edge_color = repmat(DEFAULT_EDGE_COLOR, length(elec.names), 1);
    end
    
    view(viewpoint);
    
    if(float_to_camera > 0)
        cam_vector = campos - camtarget;
        cam_vector_length = sqrt(sum(cam_vector.^2));
        elec.coords = elec.coords + repmat(float_to_camera .* cam_vector ./ cam_vector_length, size(elec.coords, 1), 1);
        if(exist('extra_float_to_camera', 'var') && ~isempty(extra_float_to_camera))
            elec.coords(extra_float_to_camera,:) = elec.coords(extra_float_to_camera,:) ...
                + 2 * repmat(cam_vector ./ cam_vector_length, size(elec.coords(extra_float_to_camera,:), 1), 1);
        end
    end
    elec.normals = compute_elec_normals(elec.coords, viewpoint, orient_to_camera);
    
    for i=1:size(elec.coords, 1)
        plot_electrode(elec.coords(i,:), elec.normals(i,:), DEFAULT_RADIUS * elec.radius(i), elec.color(i,:), elec.edge_color(i,:));
        
        if(showing_labels || (showing_signif_labels && elec.signif(i)))
            %             text('Position', [elec.coords(i, 1)+3 elec.coords(i, 2)+3 elec.coords(i, 3)], 'String', elec.names{i}, ...
            %                 'Color', round(1 - get(gca, 'Color')), 'FontSize', 8, 'Tag', 'electrode');
            text('Position', [elec.coords(i,:) + 4*camup], 'String', elec.names{i}, ...
                'Color', round(1 - get(gca, 'Color')), 'FontSize', 8, 'Tag', 'electrode');
        end
    end
    
    elec_view = @update_view_and_elecs;
    
    function update_view_and_elecs(viewpoint)
        %set(0, 'CurrentFigure', brain_figh);
        delete(findobj(gcf, 'Tag', 'electrode'));
        view(viewpoint);
        plot_electrodes(elec, viewpoint, orient_to_camera, float_to_camera, showing_labels, showing_signif_labels, extra_float_to_camera);
    end
end

function plot_electrode(center, normal, radius, elec_color, edge_color)
    persistent theta costheta sintheta;
    if(isempty(theta))
        theta = 0:0.2:2*pi;
        costheta = cos(theta);
        sintheta = sin(theta);
    end
    
    v = null(normal);
    points = repmat(center', 1, size(theta,2)) + radius * (v(:,1) * costheta + v(:,2) * sintheta);
    patch(points(1,:), points(2,:), points(3,:), elec_color, 'FaceLighting', 'none', 'EdgeColor', edge_color, 'Tag', 'electrode');
end

function normals = compute_elec_normals(coords, viewpoint, orient_to_camera)
    NEIGHB_DIST = 5;
    viewpoint_radians = pi / 180 * viewpoint;
    viewpoint_radians(1) = viewpoint_radians(1) + pi/2; %correct for different interpretations of theta betw view() and sph2cart()

    if(orient_to_camera)
        [normal(1) normal(2) normal(3)] = sph2cart(viewpoint_radians(1), viewpoint_radians(2), 1);
        normals = repmat(normal, size(coords, 1), 1);
    else
        normals = zeros(size(coords, 1), 3);
        [vertices vertex_normals] = get_fig_vertex_info();
        
        for i=1:size(coords, 1)
            wh_local_vertices = find(abs(vertices(:,1) - coords(i,1)) < NEIGHB_DIST ...
                & abs(vertices(:,2) - coords(i,2)) < NEIGHB_DIST ...
                & abs(vertices(:,3) - coords(i,3)) < NEIGHB_DIST);
            
            if(isempty(wh_local_vertices))
                warning(sprintf('Unable to find local vertices within %d of %d %d %d\n', NEIGHB_DIST, coords(i,:)));
                [normals(i,1) normals(i,2) normals(i,3)] = sph2cart(viewpoint_radians(1), viewpoint_radians(2), 1);
            else
                [theta phi] = cart2sph(vertex_normals(wh_local_vertices,1), vertex_normals(wh_local_vertices,2), vertex_normals(wh_local_vertices,3));
                [normals(i,1) normals(i,2) normals(i,3)] = sph2cart(mean(theta), mean(phi), 1);
            end
        end
    end
end
