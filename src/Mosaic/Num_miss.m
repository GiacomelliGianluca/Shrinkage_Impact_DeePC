function [miss, labels_id, labels_k] = Num_miss(str_cl_type, data_uy, conds, Centr, labels_k_p, inst, N_ini, f_inv)
% Num_miss computes the number of missclassification performed by kmeans
% concerning the future and past horizons by the ideal and k-means labels
% of its points

% - str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - data_uy: data to identify (past + future horizon) [input output]
% - conds: conditions on which we define the PWA regions (struct, where .y: output)
% - Centr: centroids of the clusters
% - inst: execution instant
% - labels_k_p: past labels of the past and future horizons points by kmeans
% - N_ini: Number of initialization steps
% - init: Used for first execution 
% - f_inv: flag for inverting the labels

% - miss: number of missclassification
% - labels_id: ideal labels of the past and future horizons points
% - labels_k: labels of the past and future horizons points by kmeans


% Extracting input data
data_u = data_uy(:,1);
% Extracting output data
data_y = data_uy(:,2);

% Extracting N_h
N_h = length(data_u);

 % Ideal labelling (1 if data_y<0, 2 if data_y>=0)
% Extracting the switching conditions 
cond_y = conds.y;
% label_id: N_h vector with the ideal labels
labels_id = NaN*ones(N_h,1);
for i = 1 : N_h
    if data_y(i) >= cond_y
        labels_id(i) = 2;
    else 
        labels_id(i) = 1;
    end
end

if strcmp(str_cl_type,'cl') % k-means case
    if inst == 1 % Initilization 
        % labels_k: N_h vector with the k-means labels
        labels_k = labels_k_p;
        for i = 1 : (length(data_uy(:,2))-N_ini)
            % Extracting data
            Data_cl = [data_u(i:i+N_ini)' data_y(i:i+N_ini-1)'];    
            % Clustering
            [~,labels_k(N_ini+i)] = pdist2(Centr,Data_cl,'euclidean','Smallest',1);
        end

        if f_inv % Switching the labels
            for i = 1 : length(labels_k)
                if labels_k(i) == 1
                    labels_k(i) = 2;
                else 
                    labels_k(i) = 1;
                end
            end
        end
    
    else % Normal case (first in, Last Out of points)
        % Extracting the last data
        Data_last = [data_u(end-N_ini:end)' data_y(end-N_ini:end-1)'];
        % Clustering the last data
        [~, label_last] = pdist2(Centr,Data_last,'euclidean','Smallest',1);
        if f_inv % Switching the last label
           if label_last == 1
              label_last = 2;
           else
                 label_last = 1;
           end
        end
        % Update of the labels
        labels_k = [labels_k_p(2:end,1);    label_last];
    end
    
    miss = sum(labels_id ~= labels_k);
else % Ideal case
    labels_k = labels_k_p;
    miss = 0;
end

end

