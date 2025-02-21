function fig_handle = Cluster_data_plot_manual(Data, label, title_string, data1_string, data2_string, fig_index)
% Cluster_data_plot realizes the 2d plot of the clustered data obtained from a
% clustering technique (e.g. kmeans)

% Data: vector of data where each column corresponds to a data type
% label: contains the label of the correspondant cluster of each data pair
% n_m: number of clusters
% title string: title of the plot
% data1_string: string x-axis (corresponds with data first column)
% data2_string: string x-axis (corresponds with data second column)
% fig_index: index of the current figure (if no figure, 0)

if fig_index > 0
    fig_handle = figure(fig_index);
    hold on, grid on,
    xlabel(data1_string,'Interpreter','latex');
    ylabel(data2_string,'Interpreter','latex');
    title(title_string, 'Interpreter','latex', 'FontWeight', 'bold');

    % col = hsv(2*n_m);
    % 
    % for i = 1 : n_m
    %     plot(Data(label==i,1),Data(label==i,2),'DisplayName', ['Cluster Mode ', num2str(i)], 'Color', col(i,:), 'Marker', '.','MarkerSize',12,'LineStyle','None')
    % end

    plot(Data(label==1,1),Data(label==1,2),'DisplayName', ['Cluster Mode ', num2str(1)], 'Color', [0.10196078431372549 0.5019607843137255 0.7333333333333333], 'Marker', '.','MarkerSize',12,'LineStyle','None') % Mode 1
    plot(Data(label==2,1),Data(label==2,2),'DisplayName', ['Cluster Mode ', num2str(2)], 'Color', [0.6274509803921569 0 0], 'Marker', '.','MarkerSize',12,'LineStyle','None') % Mode 2
    
    legend show;
    legend('Location','best');
    
else
    fig_handle = NaN;
end


end

