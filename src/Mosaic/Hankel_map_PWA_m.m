function fig_H_map = Hankel_map_PWA_m(str_cl_type, M, H_label, plt_opt, fig_index, string_matr_type)
% Function that plots the data stored in the Hankel matrix H or Mosaic matrix
% depending on their value.

% str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% M: number of modes
% - H_label: Hankel matrix of the same size of H, containing the labels
% - plt_opt: plot options

% Allocating the figure (columns on x-axis, raws on y-axis)
fig_H_map = figure(fig_index);
grid on, hold on,
str = 'Column index';
x_title = xlabel(str,'Interpreter','latex');
str_y = 'Raw index';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title(['Data Structure - ', string_matr_type, ' matrix'],'Interpreter','latex', 'FontWeight', 'bold');
xlim([1 length(H_label(1,:))]);
xticks([1 200 400 600 800]);
ylim([1 44]);
yticks([1 10 20 30 40])
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';
set(h_title, 'FontSize', plt_opt.title_font_dim);
set(y_title, 'FontSize', plt_opt.y_font_dim);
set(x_title, 'FontSize', plt_opt.x_font_dim);
fig_H_map.Units               = plt_opt.plot_unit;
fig_H_map.Position(3)         = plt_opt.plot_dim_x;
fig_H_map.Position(4)         = plt_opt.plot_dim_y;
set(fig_H_map.Children, ...
    'FontName',     plt_opt.fontname, ...
    'FontSize',     plt_opt.fontsize-2);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_H_map.PaperPositionMode   = 'auto';

% Specify the size
sz = 4;

% Defining the modes
modes = 1:M;

% Identifying the coordinates in the matrix 
[raws, col] = find(H_label);

% Extraction of the labels
label = H_label(sub2ind(size(H_label), raws, col));



if strcmp(str_cl_type,'id')
    
    scatter(col(label == 1), raws(label == 1), sz, [0.10196078431372549 0.5019607843137255 0.7333333333333333],'filled'); % data mode_1
    scatter(col(label == 2), raws(label == 2), sz, [0.6274509803921569 0 0],'filled');  % data mode_2
 
    legend(['Mode ', num2str(modes(1))], ['Mode ', num2str(modes(2))], 'Location','best');
end

if strcmp(str_cl_type,'cl')

    scatter(col(label == 1), raws(label == 1), sz, [0.10196078431372549 0.5019607843137255 0.7333333333333333],'filled'); % data mode_1
    scatter(col(label == 2), raws(label == 2), sz, [0.6274509803921569 0 0],'filled');  % data mode_2

    legend(['Mode ', num2str(modes(1))], ['Mode ', num2str(modes(2))], 'Location','best');

end



end

