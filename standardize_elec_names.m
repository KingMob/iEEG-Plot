function label = standardize_elec_names(label)
    convert_to_char = false;
    if(ischar(label))
        label = cellstr(label);
        convert_to_char = true;
    end
    
    %tic;
    for i=1:length(label(:))
        elec_parts = textscan(label{i}, '%[^0123456789]%d');
        label{i} = sprintf('%s%02d', elec_parts{1}{1}, elec_parts{2});
    end
    %toc
    
    if(convert_to_char)
        label = char(label);
    end
end

% function newname = standardize_elec_name(oldname)
%     elec_parts = textscan(oldname, '%[^0123456789]%d');
%     newname = sprintf('%s%02d', elec_parts{1}{1}, elec_parts{2});
% end
