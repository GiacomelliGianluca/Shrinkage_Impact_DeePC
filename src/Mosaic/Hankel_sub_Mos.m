function [sub_ind_H, fig_H_sub] = Hankel_sub_Mos(mode_M, N_h, n_x, n_u, n, fig_H, plt_opt, fig_index)
% Function that identifies the submatrixes associated to the different
% modes in the Hankel matrix H 

% - sub_ind_data: matrix of the indexes of the submatrices associated with the data
% (where the data considered are the ones of the first column and the
% last raw (!! hence, with mosaic not all data are covered and there is not
% an allignment with the database indexes!!)). See below for detail
% - sub_ind_H: matrix of the indexes associated with the position (raw and column)
% in the Hankel matrix of the submatrices. See below for detail
% - sub_num: number of subamtrices
% - fig_H_sub: handle for the submatrix figure

% - mode_M: matrix of the same size of H_M, containing the
% mode (centroid) of the partitioned data
% - N_h: Hankel matrix depth
% - n_x: state degree
% - n_u: input degree
% - n: number of Hankel matrixes. It corresponds to the number of clusters
% - fig_H: Hankel (or Mosaic) figure handle

% Indexes initialization
% -Column 1: starting index -Column 2: ending index
% -Column 3: belonging mode -Column 4: SP_id_cond (1: ok, 0: not ok)
% Each raw describe a differen submatrices
sub_ind_H = [];

% sub_num: counter for the number of submatrices. (it will be equal to the
% selected number of clusters n)
sub_num = 1;

% Initialization
sub_ind_H(sub_num, 1) = 1;
sub_ind_H(sub_num, 3) = mode_M(1,1);

for i = 2 : length(mode_M(1,:))
    
    if mode_M(1,i) ~= mode_M(1,i-1)
        % Saving last index of the previous mode
        sub_ind_H(sub_num, 2) = i-1;
        % Changing mode and saving the initial index and its value
        sub_num = sub_num + 1;
        sub_ind_H(sub_num, 1) = i;
        sub_ind_H(sub_num, 3) = mode_M(1,i);
    end    

end

% Check
if sub_num ~= n
    keyboard; %error
end

% Finalization
sub_ind_H(n, 2) = length(mode_M(1,:));


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
    allChildren = get(gca, 'Children');                % list of all objects on axes
    displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
    delete(allChildren(1:length(displayNames)/2));
    grai = gray(sub_num);
    for j = 1 : sub_num
       xline(sub_ind_H(j,1), 'DisplayName', ['Submatrix ', num2str(j),' - begin'], 'Color', grai(j,:), LineWidth=2)
       xline(sub_ind_H(j,2), 'DisplayName', ['Submatrix ', num2str(j),' - end'], 'Color', grai(j,:), LineWidth=2)
    end
    legend('NumColumns',7)
end

%% Checking the SP id conditions for each submatrix
% Subspace predictor id conditions
SP_id_cond = (N_h + n_x)*(n_u + 1) - 1;

% For each of the submatrices, we verified if we have enough data (if the persistency of excitation is satisfied, then we have enough data)
for k = 1 : sub_num
    sub_ind_H(k,4) = ((length(mode_M(:,sub_ind_H(k, 1))) + length(mode_M(end,sub_ind_H(k, 1) + 1:sub_ind_H(k, 2)))) >= SP_id_cond);
end

end

