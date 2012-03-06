function [brain_figh brain_view] = plot_brain(surf_brain, viewpoint, varargin)
    FSL_DISPLACEMENT = [-1 -17 19];
    DEFAULT_BRAIN_COLOR = [.7 .7 .7];
    
    use_fsl_displacement = false;
    brain_color = DEFAULT_BRAIN_COLOR;
    
    parse_inputs(varargin);
    
    if(~exist('viewpoint', 'var') || isempty(viewpoint)), viewpoint = [180 90]; end
    if(ischar(surf_brain)), surf_brain = cellstr(surf_brain); end
    
    if(iscellstr(surf_brain))
        surf_files = surf_brain;
        clear surf_brain;
        
        for i=1:length(surf_files)
            [d f e] = fileparts(surf_files{i});
            if(strcmp(e, '.mat'))
                surf_brain(i) = load(surf_files{i});
            else
                % See if it's a Freesurfer surf file
                try
                    [s.vertices s.faces] = read_surf(surf_files{i});
                    s.vertices = s.vertices + repmat(FSL_DISPLACEMENT, size(s.vertices, 1), 1); %correction for FSL weirdness
                    s.faces = s.faces + 1; % switch from 0-based to 1-based indexing
                    surf_brain(i) = s; clear s;
                catch exception
                    error('Unknown surface filetype in: %s\nError: %s\n', surf_files{i}, exception.message());
                    fclose('all');
                end
            end
        end
    end
    
%    delete(findobj(gcf, 'type', 'light'));
    
    for i=1:length(surf_brain)
        if(isfield(surf_brain(i), 'vertices'))
            vertices = surf_brain(i).vertices;
            faces = surf_brain(i).faces;
        elseif(isfield(surf_brain(i), 'bnd'))
            vertices = surf_brain(i).bnd.pnt;
            faces = surf_brain(i).bnd.tri;
        elseif(isfield(surf_brain(i), 'coords'))
            vertices = surf_brain(i).coords;
            faces = surf_brain(i).faces;
        end
        trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
            'FaceVertexCData', repmat(brain_color(:)', [size(vertices, 1) 1]), ...
            'FaceColor', 'interp', ...
            'FaceVertexAlphaData', .2 * ones(size(vertices, 1), 1), ...
            'Tag', 'brain');
        hold on;
    end
    
    
    % Set viewing params
    shading interp;
    lighting gouraud;
    material dull;
    axis off
    axis equal;
    camproj('orthographic');
    %hold on;
    
    update_view_and_light(viewpoint);
%    view(viewpoint);
 %   camlight('headlight', 'infinite');
    set(gcf, 'color', 'white');
    brain_figh = gcf();
    brain_view = @update_view_and_light;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function parse_inputs(args)
        if(~exist('args', 'var') || isempty(args)), return; end
        for i=1:length(args)
            if(ischar(args{i}))
                switch(lower(args{i}))
                    case 'use_fsl_displacement'
                        use_fsl_displacement = args{i+1};
                    case 'brain_color'
                        brain_color = args{i+1};
                        assert(numel(brain_color) == 3);
                    case 'translucency'
                        brain_color = args{i+1};
                        assert(numel(brain_color) == 3);
                end
            end
        end
    end
    
    function update_view_and_light(viewpoint)
        view(viewpoint);
        %delete(findobj(gcf, 'type', 'light'));
        delete(findobj(gca, 'type', 'light'));
        camlight('headlight', 'infinite');
    end
end