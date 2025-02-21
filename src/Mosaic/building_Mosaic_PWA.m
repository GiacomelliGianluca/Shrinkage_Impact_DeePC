function [H_M, index_M, label_M, mode_M, n_new, label_del] = building_Mosaic_PWA(Data, N_h, n, label, centr, d, n_x, flag_PE)
% building_Mosaic is a script to construct a mosaic Hankel matrix formed
% by n Hankel matrixes with depth N_h; if the data are not enough or
% Persistently Exciting, they are discarded.

% - N_h: past + prediction horizon (depth Hankel matrix)
% - n: number of Hankel matrixes
% - label: label of data
% - d: data dimensionality (corresponing to the number of elements for an
% istant t (e.g. number of inputs)
% - centr: matrix containing the centroids of the data
% - n_x: state degree

% - H_M: Mosaic Hankel matrix
% - index_M: matrix of the same size of H_M, containing the
% indexes of the trajectories in H_M
% - label_M: matrix of the same size of H_M, containing the
% labels of the partitioned data
% - mode_M: matrix of the same size of H_M, containing the
% mode (centroid) of the partitioned data
% n_new: new number of modes. It can be decreased if a cluster (and so the
% correspondant mosaic) has not the minimum number of data
% label_del: indeces of the delated modes

% Allocating the mosaic Hankel
H_M = [];
% Allocating the matrix with the indexes
index_M = [];
% Allocating the matrix with the labels
label_M = [];
% Allocating the matrix with the modes
mode_M = [];
% Initializing the number of modes
n_new = n;
% Vector including the index of the delated modes
label_del = [];

for i = 1 : n
    % Data_i: data having label_i
    Data_i = Data(label==i);
    % Chek #data cluster_i
    if length(Data_i) < N_h
        f = 1;
        keyboard; % We do not have enough data to create the submatrix
    else
        % index_i: indexes of the data having label_i
        index_i = find(label==i);
        % H_i: ith-Hankel matrix
        [H_i, index_i, f] = building_Hankel_w_index(Data_i, index_i, N_h, d, n_x, flag_PE);
        % Allocating the H_i Hankel into the mosaic
        H_M = [ H_M H_i];
        % Allocating the matrix of the indexes 
        index_M = [ index_M index_i];
        % Allocating the matrix of the labels
        label_M = [ label_M i*ones(size(H_i))];
        % Allocating the matrix of the centroids
        mode_M = [ mode_M centr(i)*ones(size(H_i))];
    end

    % Checking if the modes (and associated data) was delated 
    if f
        % Updating the number of modes
        n_new = n_new - f;
        label_del = [label_del   i];
    end
end

end

