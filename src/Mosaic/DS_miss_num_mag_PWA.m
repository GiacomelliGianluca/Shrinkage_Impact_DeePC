function [miss_num, miss_mag_avg_u, miss_mag_avg_y, fig_miss_num, fig_miss_mag_avg] = DS_miss_num_mag_PWA(str_cl_type, DS_label_miss, DS_u, DS_y, plt_opt, fig_index_1, fig_index_2)
% Function that plots number of misclassification for the column of a data
% structure and their avarage magnitude 

% - str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - DS_label_miss: matrix reflecting the data structure 
%                   with labeled data (misclassification is 0)
% - DS_u: input data structure
% - DS_y: output data structure
% - plt_opt: plot options

% - miss_num: number of misclassification
% - miss_mag_avg_u: input average magnitude for each column of the data
% structure
% - miss_mag_avg_y: output average magnitude for each column of the data
% structure

if ~strcmp(str_cl_type, 'id')
    miss_num = sum(DS_label_miss == 0,1);
    
    % Initializing 
    miss_mag_avg_u = NaN*zeros(length(DS_u(1,:)),1);
    miss_mag_avg_y = NaN*zeros(length(DS_u(1,:)),1);
    % Initializing correct data
    corr_mag_avg_u = NaN*zeros(length(DS_u(1,:)),1);
    corr_mag_avg_y = NaN*zeros(length(DS_u(1,:)),1);
    
    
    for i = 1 : length(DS_u(1,:))
        miss_mag_avg_u(i) = mean(abs(DS_u((DS_label_miss(:,i) == 0),i)), 1);
        miss_mag_avg_y(i) = mean(abs(DS_y((DS_label_miss(:,i) == 0),i)), 1);
        corr_mag_avg_u(i) = mean(abs(DS_u((DS_label_miss(:,i) ~= 0),i)), 1);
        corr_mag_avg_y(i) = mean(abs(DS_y((DS_label_miss(:,i) ~= 0),i)), 1);
    end
    
    
    if fig_index_1~=0
        % Allocating the figure (columns on x-axis)
        fig_miss_num = figure(fig_index_1);
        grid on, hold on,
        plot(1:length(DS_u(1,:)), miss_num, 'Color',[0.388 0.584 0])
        xlim([1 length(DS_y(1,:))]);
        str = 'Column index';
        x_title = xlabel(str,'Interpreter','latex');
        str_y = 'Number';
        y_title = ylabel(str_y,'Interpreter','latex');
        h_title = title('Number of misclassification per column','Interpreter','latex', 'FontWeight', 'bold');
        set(h_title, 'FontSize', plt_opt.title_font_dim);
        set(y_title, 'FontSize', plt_opt.y_font_dim);
        set(x_title, 'FontSize', plt_opt.x_font_dim);
        fig_miss_num.Units               = plt_opt.plot_unit;
        fig_miss_num.Position(3)         = plt_opt.plot_dim_x;
        fig_miss_num.Position(4)         = plt_opt.plot_dim_y;
        set(fig_miss_num.Children, ...
            'FontName',     plt_opt.fontname, ...
            'FontSize',     plt_opt.fontsize-2);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
        fig_miss_num.PaperPositionMode   = 'auto';
    end
    
    if fig_index_2~=0
        % Allocating the figure (columns on x-axis)
        fig_miss_mag_avg = figure(fig_index_2);
        grid on, hold on,
        plot(1:length(DS_u(1,:)), miss_mag_avg_u)
        plot(1:length(DS_y(1,:)), miss_mag_avg_y)
        plot(1:length(DS_u(1,:)), corr_mag_avg_u, 'Color','k')
        plot(1:length(DS_y(1,:)), corr_mag_avg_y, 'Color', [0.5000    0.5000    0.5000])
        str = 'Column index';
        x_title = xlabel(str,'Interpreter','latex');
        xlim([1 length(DS_y(1,:))]);
        str_y = 'Magnitude';
        y_title = ylabel(str_y,'Interpreter','latex');
        h_title = title('Average magnitude of misclassification','Interpreter','latex', 'FontWeight', 'bold');
        set(h_title, 'FontSize', plt_opt.title_font_dim);
        set(y_title, 'FontSize', plt_opt.y_font_dim);
        set(x_title, 'FontSize', plt_opt.x_font_dim);
        fig_miss_mag_avg.Units               = plt_opt.plot_unit;
        fig_miss_mag_avg.Position(3)         = plt_opt.plot_dim_x;
        fig_miss_mag_avg.Position(4)         = plt_opt.plot_dim_y;
        set(fig_miss_mag_avg.Children, ...
            'FontName',     plt_opt.fontname, ...
            'FontSize',     plt_opt.fontsize-2);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
        fig_miss_mag_avg.PaperPositionMode   = 'auto';
        legend('Miss input','Miss output', 'Correct input','Correct output', 'Location','best')
    end
else
    miss_num = NaN;
    miss_mag_avg_u = NaN;
    miss_mag_avg_y = NaN;
    fig_miss_num = NaN;
    fig_miss_mag_avg = NaN;
end





end

