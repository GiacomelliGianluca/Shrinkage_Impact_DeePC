function [sub_ind_data_miss, sub_ind_H_miss, sub_num_miss, fig_H_sub_miss] = Hankel_sub_miss_PWA(str_cl_type, Matr_label_miss, sub_ind_H, sub_num, H_u, N_h, n_x, n_u, fig_H_miss,  plt_opt, fig_index)
% Function that identifies the submatrixes associated to the different
% modes in the Hankel matrix H 

% - str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - sub_ind_data_miss: matrix of the indeces of the submatrices associated with the data
% (where the data considered are the ones of the first column and the
% last raw (!! hence, with mosaic not all data are covered and there is not
% an alignment with the database indexes!!))
% - sub_ind_H_miss: matrix of the indexes associated with the position (raw and column)
% in the Hankel matrix of the submatrices. See below for detail
% - sub_num_miss: number of subamtrices without misclassification
% - fig_H_sub: handle for the submatrix figure

% - Matr_label_miss: matrix reflecting the Hankel structure 
%                   with labeled data (misclassification is 0)
% - Matr_index: matrix reflecting the Hankel structure containing the indeces
%            of the data
% - sub_ind_H: matrix of the indexes associated with the position (raw and column)
% in the Hankel matrix of the submatrices. See below for detail
% - sub_num: number of subamtrices
% - H: Hankel matrix
% - n_x: state degree
% - n_u: input degree
% - fig_H: Hankel (or Mosaic) figure handle

if ~strcmp(str_cl_type, 'id')
    % Subspace predictor id conditions
    SP_id_cond = (N_h + n_x)*(n_u + 1) - 1;
    
    % Retrieving the data vector
    data_u = [H_u(:,1)'   H_u(end,2:end)];
    
    % Indexes initialization
    % -Column 1: starting index -Column 2: ending index -Column 3: belonging mode
    % Each raw describe a differen submatrices
    sub_ind_data_miss = [];
    
    % f_b_submatr: flag identifying if a submatrix is being identify
    % 0: searching for beginning index
    % 1: searching for final index (mode 1)
    % 2: searching for final index (mode 2)
    % 4: searching for final index (mode 3)
    % 6: searching for final index (mode 4)
    % ...
    % n: searching for final index (mode n)
    f_b_submatr = 0;
    
    % sub_num_miss: counter for the number of submatrices
    sub_num_miss = 0;
    
    % if a matrix is found, check SP_id_cond and PE
    % Identification of the submatrices associated to the mode 1
    for i = 1 : sub_num
       % Retrieving the label vector
       labels = [Matr_label_miss(:,sub_ind_H(i,1))'   Matr_label_miss(end,sub_ind_H(i,1)+1:sub_ind_H(i,2))];
       for j = 1 : (length(labels) - SP_id_cond + 1 )
           if (sum(labels(j : (j - 1) + SP_id_cond) > 0) == SP_id_cond) || (f_b_submatr == 1)
               if f_b_submatr == 0
                   sub_num_miss = sub_num_miss + 1;
                   f_b_submatr = 1;
                   sub_ind_data_miss(sub_num_miss, 1) = sub_ind_H(i,1) + j - 1;
                   sub_ind_data_miss(sub_num_miss, 3) = sub_ind_H(i,3);
               end
               %
               if sum(labels(j : (j - 1) + SP_id_cond) > 0) ~= SP_id_cond
                   f_b_submatr = 0;
                   sub_ind_data_miss(sub_num_miss, 2) = sub_ind_H(i,1) - 1 + (j - 2) + SP_id_cond;
               end
               %
               if j == (length(labels) - SP_id_cond + 1) && (f_b_submatr == 1)
                   f_b_submatr = 0;
                   sub_ind_data_miss(sub_num_miss, 2) = sub_ind_H(i,2) + N_h - 1;
               end
           end           
       end
    end
    
    % Checking persistency of excitation of the submatrices
    for k = 1 : sub_num_miss
       PE_check(data_u(sub_ind_data_miss(k,1):sub_ind_data_miss(k,2)), (n_u* (N_h + n_x)));
    end
    
    % Plotting (the final indexes)
    if sub_num_miss > 0
        % Allocating the figure (columns on x-axis, raws on y-axis)
        fig_H_sub_miss = figure(fig_index); hold on
        CloneFig(fig_H_miss, fig_H_sub_miss);
        ax = gca;
        ax.YDir = 'reverse';
        ax.XAxisLocation = 'top';
        fig_H_sub_miss.Units               = plt_opt.plot_unit;
        fig_H_sub_miss.Position(3)         = plt_opt.plot_dim_x;
        fig_H_sub_miss.Position(4)         = plt_opt.plot_dim_y;
        set(fig_H_sub_miss.Children, ...
            'FontName',     plt_opt.fontname, ...
            'FontSize',     plt_opt.fontsize-2);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
        fig_H_sub_miss.PaperPositionMode   = 'auto';
        % Adjusting the legend
        allChildren = get(gca, 'Children');                % list of all objects on axes
        displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
        delete(allChildren(1:length(displayNames)/2));
        grai = gray(sub_num_miss);
        for j = 1 : sub_num_miss
           xline(sub_ind_data_miss(j,1), 'DisplayName', ['Submatrix ', num2str(j),' - begin'], 'Color', grai(j,:), LineWidth=2)
           xline(sub_ind_data_miss(j,2) - N_h + 1, 'DisplayName', ['Submatrix ', num2str(j),' - end'], 'Color', grai(j,:), LineWidth=2)
        end
        legend('NumColumns',(2*length(unique(sub_ind_H(:,3)))-1))
    
        % Identify Hankel indexes 
        sub_ind_H_miss(:,1) = sub_ind_data_miss(:, 1);
        sub_ind_H_miss(:,3) = sub_ind_data_miss(:, 3);
        for k = 1 : sub_num_miss
            sub_ind_H_miss(k, 2) = sub_ind_data_miss(k, 2) - N_h + 1;
        end
    else 
        fig_H_sub_miss = NaN;
        sub_ind_H_miss = sub_ind_data_miss;
    end
else
    sub_ind_data_miss = NaN;
    sub_ind_H_miss = NaN;
    sub_num_miss = NaN;
    fig_H_sub_miss = NaN;
end

