function [M_label_miss, fig_H_map_miss] = Hankel_map_miss_PWA(str_cl_type, fig_H, label, label_id, f_inv, M_label, M_index, N_ini, plt_opt, fig_index, string_matr_type)
% Function that plots the data stored in the Hankel matrix H or Mosaic matrix
% depending on their value highlighting the missclassfied data. 

% - str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - fig_H: Hankel (or Mosaic) figure handle
% - label: labels assigned by the clustering technique
% - label_id: ideal labels
% - f_inv: flag for inverting the labels
% - M_label: matrix reflecting containing the labeled data structure 
% - M_index: matrix reflecting the Hankel structure containing the indeces
%            of the data


% label_miss: dataset pointing out the missclassfied data


% Specify the size
sz = 4;

if ~strcmp(str_cl_type, 'id')
    % Allocating the figure (columns on x-axis, raws on y-axis)
    fig_H_map_miss = figure(fig_index); hold on
    CloneFig(fig_H, fig_H_map_miss);
    ax = gca;
    ax.YDir = 'reverse';
    ax.XAxisLocation = 'top';
    fig_H_map_miss.Units               = plt_opt.plot_unit;
    fig_H_map_miss.Position(3)         = plt_opt.plot_dim_x;
    fig_H_map_miss.Position(4)         = plt_opt.plot_dim_y;
    set(fig_H_map_miss.Children, ...
        'FontName',     plt_opt.fontname, ...
        'FontSize',     plt_opt.fontsize-2);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
    fig_H_map_miss.PaperPositionMode   = 'auto';
    % Adjusting the legend
    allChildren = get(gca, 'Children');                
    displayNames = get(allChildren, 'DisplayName');    
    delete(allChildren(1:length(displayNames)/2));
    
    if f_inv
        for i = 1 : length(label_id)
            if label_id(i) == 1
               label_id(i) = 2;
            else 
                label_id(i) = 1;                
            end
        end
    end

    % label_miss: missclassfied data
    label_miss = label ~= label_id;
    % Identifying the coordinates in the matrix (you cannot use H, becuase it
    % has zeroes and so find does not recognize these entries)
    [raws, col] = find(M_index);
    % Extraction of the data corresponding to the coordinates
    data = M_index(sub2ind(size(M_index), raws, col));

    % M_label_miss: initialization
    M_label_miss = M_label;

    if strcmp(string_matr_type, 'Hankel')
        % f_legend: flag for the legend
        f_legend = 1;
        % label index
        label_index = 1 : length(label_id);
        for i = 1 : length(label_miss)
            if label_miss(i) == 1
                % Changing the label as misclassified 
                M_label_miss(M_index == i) = 0;
                if f_legend
                    scatter(col(data == label_index(i)), raws(data == label_index(i)), sz, [0.388 0.584 0],'filled', 'DisplayName', 'Missclassified'); 
                    f_legend = 0;
                else
                    scatter(col(data == label_index(i)), raws(data == label_index(i)), sz, [0.388 0.584 0],'filled', 'HandleVisibility','off'); 
                end

            end
        end
   
    end

    if strcmp(string_matr_type, 'Mosaic') % Same as Hankel, where the first N_ini data are not considered
        % f_legend: flag for the legend
        f_legend = 1;
        % label index
        label_index = 1 : length(label_id);
        for i = N_ini + 1 : length(label_miss)
            if label_miss(i) == 1
                % Changing the label as misclassified 
                M_label_miss(M_index == i) = 0;
                if f_legend
                    scatter(col(data == label_index(i)), raws(data == label_index(i)), sz, 'magenta','filled', 'DisplayName', 'Missclassified'); 
                    f_legend = 0;
                else
                    scatter(col(data == label_index(i)), raws(data == label_index(i)), sz, 'magenta','filled', 'HandleVisibility','off'); 
                end

            end
        end
   
    end

    
else % Ideal case
     M_label_miss = M_label;
     fig_H_map_miss = NaN;
end


end

