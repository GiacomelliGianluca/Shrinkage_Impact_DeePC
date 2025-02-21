function [sub_ind_data, sub_ind_H, sub_num, fig_H_sub] = Hankel_sub_PWA(str_cl_type, H_y, H_u, M, N_h, n_x, n_u, conds, H_label, fig_H,  plt_opt, fig_index)
% Function that identifies the submatrixes associated to the different
% modes in the Hankel matrix H 

% - sub_ind_data: matrix of the indexes of the submatrices associated with the data
% (where the data considered are the ones of the first column and the
% last raw (!! hence, with mosaic not all data are covered and there is not
% an alignment with the database indexes!!)). See below for detail
% - sub_ind_H: matrix of the indexes associated with the position (raw and column)
% in the Hankel matrix of the submatrices. See below for detail
% - sub_num: number of subamtrices
% - fig_H_sub: handle for the submatrix figure

% str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - H: Hankel matrix
% - M: number of modes
% - conds: conditions on which we define the PWA regions (struct, where .y: output)
% - n_x: state degree
% - n_u: input degree
% - H_label: Hankel matrix of the same size of H, containing the labels
% - fig_H: Hankel (or Mosaic) figure handle


% Defining the modes
modes = 1:M;

% Subspace predictor id conditions
SP_id_cond = (N_h + n_x)*(n_u + 1) - 1;

% Retrieving the data vector
data_y = [H_y(:,1)'   H_y(end,2:end)];
data_u = [H_u(:,1)'   H_u(end,2:end)];
 

% Indexes initialization
% -Column 1: starting index -Column 2: ending index -Column 3: belonging mode
% Each raw describe a differen submatrices
sub_ind_data = [];

% f_b_submatr: flag identifying if a submatrix is being identify
% 0: searching for beginning index
% 1: searching for final index (mode 1)
% 2: searching for final index (mode 2)
% 4: searching for final index (mode 3)
% 6: searching for final index (mode 4)
% ...
% n: searching for final index (mode n)
f_b_submatr = 0;

% sub_num: counter for the number of submatrices
sub_num = 0;

if strcmp(str_cl_type,'id')
    % Extracting the switching conditions 
    cond_y = conds.y;
    
    for i = 1 : (length(data_y) - SP_id_cond + 1 )
        
        % Identification of the submatrices associated to the mode 1
        if (sum(data_y(i : (i - 1) + SP_id_cond) < cond_y) == SP_id_cond) || (f_b_submatr == 1)
           if f_b_submatr == 0
               sub_num = sub_num + 1;
               f_b_submatr = 1;
               sub_ind_data(sub_num, 1) = i;
               sub_ind_data(sub_num, 3) = modes(1);
           end
           %
           if sum(data_y(i : (i - 1) + SP_id_cond) < cond_y) ~= SP_id_cond
               f_b_submatr = 0;
               sub_ind_data(sub_num, 2) = (i - 2) + SP_id_cond;
           end
           %
           if i == (length(data_y) - SP_id_cond + 1) && (f_b_submatr == 1)
               f_b_submatr = 0;
               sub_ind_data(sub_num, 2) = length(data_y);
           end
        end
    
        % Identification of the submatrices associated to the mode M
        if M > 1
            if (sum(data_y(i : (i - 1) + SP_id_cond) > cond_y) == SP_id_cond) || (f_b_submatr == M)
               if f_b_submatr == 0
                   sub_num = sub_num + 1;
                   f_b_submatr = M;
                   sub_ind_data(sub_num, 1) = i;
                   sub_ind_data(sub_num, 3) = modes(end);
               end
               %
               if sum(data_y(i : (i - 1) + SP_id_cond) > cond_y) ~= SP_id_cond
                   f_b_submatr = 0;
                   sub_ind_data(sub_num, 2) = (i - 2) + SP_id_cond;
               end
            end
        end

       % Exit case
       if i == (length(data_y) - SP_id_cond + 1) && (f_b_submatr > 0)
           f_b_submatr = 0;
           sub_ind_data(sub_num, 2) = length(data_y);
       end
    
    end
end

if strcmp(str_cl_type,'cl')
    % Retrieving the label vector
    label = [H_label(:,1)'   H_label(end,2:end)];
    for i = 1 : (length(label) - SP_id_cond + 1 )
        
        % Identification of the submatrices associated to the mode 1
        if (sum(label(i : (i - 1) + SP_id_cond) == 1) == SP_id_cond) || (f_b_submatr == 1)
           if f_b_submatr == 0
               sub_num = sub_num + 1;
               f_b_submatr = 1;
               sub_ind_data(sub_num, 1) = i;
               sub_ind_data(sub_num, 3) = modes(1);
           end
           %
           if sum(label(i : (i - 1) + SP_id_cond) ==1) ~= SP_id_cond
               f_b_submatr = 0;
               sub_ind_data(sub_num, 2) = (i - 2) + SP_id_cond;
           end
           %
           if i == (length(label) - SP_id_cond + 1) && (f_b_submatr == 1)
               f_b_submatr = 0;
               sub_ind_data(sub_num, 2) = length(label);
           end
        end
    
        % Identification of the submatrices associated to the mode M
        if M > 1
            if (sum(label(i : (i - 1) + SP_id_cond) == M) == SP_id_cond) || (f_b_submatr == M)
               if f_b_submatr == 0
                   sub_num = sub_num + 1;
                   f_b_submatr = M;
                   sub_ind_data(sub_num, 1) = i;
                   sub_ind_data(sub_num, 3) = modes(end);
               end
               %
               if sum(label(i : (i - 1) + SP_id_cond) == M) ~= SP_id_cond
                   f_b_submatr = 0;
                   sub_ind_data(sub_num, 2) = (i - 2) + SP_id_cond;
               end
            end
        end
    
       % Exit case
       if i == (length(label) - SP_id_cond + 1) && (f_b_submatr > 0)
           f_b_submatr = 0;
           sub_ind_data(sub_num, 2) = length(label);
       end
    
    end
end


% Checking persistency of excitation of the submatrices
for j = 1 : sub_num
   PE_check(data_u(sub_ind_data(j,1):sub_ind_data(j,2)), (n_u* (N_h + n_x)));
end

% Plotting (the final indexes)
if sub_num > 0
    % Allocating the figure (columns on x-axis, raws on y-axis)
    fig_H_sub = figure(fig_index); hold on
    CloneFig(fig_H, fig_H_sub);
    ax = gca;
    ax.YDir = 'reverse';
    ax.XAxisLocation = 'top';
    fig_H_sub.Units               = plt_opt.plot_unit;
    fig_H_sub.Position(3)         = plt_opt.plot_dim_x;
    fig_H_sub.Position(4)         = plt_opt.plot_dim_y;
    set(fig_H_sub.Children, ...
        'FontName',     plt_opt.fontname, ...
        'FontSize',     plt_opt.fontsize-2);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
    fig_H_sub.PaperPositionMode   = 'auto';
    % Adjusting the legend
    allChildren = get(gca, 'Children');                
    displayNames = get(allChildren, 'DisplayName');    
    delete(allChildren(1:length(displayNames)/2));
    grai = gray(sub_num);
    for j = 1 : sub_num
       xline(sub_ind_data(j,1), 'DisplayName', ['Submatrix ', num2str(j),' - begin'], 'Color', grai(j,:), LineWidth=2)
       xline(sub_ind_data(j,2) - N_h + 1, 'DisplayName', ['Submatrix ', num2str(j),' - end'], 'Color', grai(j,:), LineWidth=2)
    end
    legend('NumColumns',(2*M)-1)

    % Identify Hankel indexes 
    sub_ind_H(:,1) = sub_ind_data(:, 1);
    sub_ind_H(:,3) = sub_ind_data(:, 3);
    for k = 1 : sub_num
        sub_ind_H(k, 2) = sub_ind_data(k, 2) - N_h + 1;
    end
else 
    fig_H_sub = NaN;
    sub_ind_H = sub_ind_data;
end

end

