function fig_BPI = plotting_BPI_single_h(BPI_cell, plt_opt, fig_index)

% - BPI_cell: cell strcuture containing in each row the bhevaioral
% performance index (BPI) of a mode. The element has then as rows the time instants and as columns
% the hyperparameter value
% - lambda_g_x: chosen lambda_g
% - plt_opt: plot options

% - fig_BPI: plot of the BPI of each mode in the different control instants

str_x = 'Step';
str_y = '$BPI_{i,t}$';


% Plotting weights in the different time instants
if fig_index ~=0
   fig_BPI = figure(fig_index);
   hold on, grid on, 
   color = [0.10196078431372549 0.5019607843137255 0.7333333333333333;
            0.6274509803921569 0 0;
            0 0 0];
   plt = NaN*ones(length(BPI_cell(:,1)),1);
   for i = 1 : length(BPI_cell(:,1))
      if i == length(BPI_cell(:,1)) % BPI none modes case
        plt(i)=plot(0:1:length(BPI_cell{i,1}(:,1))-1, BPI_cell{i,1}(:,1), 'Color', color(i,:), 'linewidth',plt_opt.thick,'DisplayName', 'Mode None');
      else % Normal case
        plt(i)=plot(0:1:length(BPI_cell{i,1}(:,1))-1, BPI_cell{i,1}(:,1), 'Color', color(i,:), 'linewidth',plt_opt.thick,'DisplayName', ['Mode ', num2str(i)]);
      end
   end
   x_title_b = xlabel(str_x,'Interpreter','latex');
   y_title_b = ylabel(str_y,'Interpreter','latex');
   str_title = 'BPI';
   h_title_b = title(str_title, 'Interpreter','latex', 'FontWeight', 'bold');
   yline(1, "--", 'Color', [0 0 0], 'DisplayName', '');
   xlim([0 49]);
   legend (plt)
   legend( 'Location','best','Interpreter','latex');
   set(h_title_b, 'FontSize', plt_opt.title_font_dim);
   set(y_title_b, 'FontSize', plt_opt.y_font_dim);
   set(x_title_b, 'FontSize', plt_opt.x_font_dim);
   fig_BPI.Units               = plt_opt.plot_unit;
   fig_BPI.Position(3)         = plt_opt.plot_dim_x;
   fig_BPI.Position(4)         = plt_opt.plot_dim_y;
   set(fig_BPI.Children, ...
        'FontName',     plt_opt.fontname, ...
        'FontSize',     plt_opt.fontsize-2);
   set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
   fig_BPI.PaperPositionMode   = 'auto';
end

end

