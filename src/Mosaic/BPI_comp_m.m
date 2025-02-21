function BPI_singles = BPI_comp_m(str_cl_type, g, ind_H, n_x, n_u, n_m, labels, inst, N_ini, f_inv)
% BPI_comp_mos computes the Behavioural Performance Index (BPI) for an instant 

% HANKEL VERSION

% - str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - g: selector
% - ind_H: indexes of the submatrices -Column 1: starting index -Column 2: ending index
% - n_x: system state
% - n_u: numeber of inputs
% - n_m: number of modes
% - inst: execution instant
% - labels: labels of the past and future horizons points
% - N_ini: Number of initialization steps
% - f_inv: flag for inverting the labels



% thr: sensibility of the zero norm
thr = 1e-2;

if f_inv && (strcmp(str_cl_type,'cl')) % Switching the labels to be coherent with cl (previosuly they were coherent with the ideal one, but not with the g values)
    for i = 1 : length(labels)
        if labels(i) == 1
           labels(i) = 2;
        else
            labels(i) = 1;
        end
    end
end

if strcmp(str_cl_type,'id') || (strcmp(str_cl_type,'cl') && (inst >= N_ini + 1))
    % Allocating BPI_singles
    BPI_singles = NaN*ones(n_m,1);
    for i = 1 : n_m
        % mode_pos_index: Index pointing the position of mode_i in ind_H
        mode_pos_index = find(ind_H(:,3) == i);
        % N_h_i: number of elements associated with the mode_i in N_h
        N_h_i = sum(labels == i);
        % n_g: number of not null elements in g associated with the mode i
        n_g = sum( abs(g(ind_H(mode_pos_index,1) : ind_H(mode_pos_index,2))) >= thr);
        BPI_singles(mode_pos_index) = n_g/(n_u*N_h_i+n_x);
    end
else 
        BPI_singles = NaN*ones(n_m,1);
end

end

