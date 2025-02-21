function BPI_singles = BPI_comp_h(str_cl_type, g, ind_H, num_sub , n_x, n_u, n_m, labels, inst, N_ini, f_inv)
% BPI_comp_mos computes the Behavioural Performance Index (BPI) for an instant (avarage of modes in BPI_avg and single modes in BPI_singles)

% iHANKEL VERSION

% - str_cl_type: clustering type ('id': ideal case, 'cl': kmeans)
% - g: selector
% - ind_H: indexes of the submatrices -Column 1: starting index -Column 2: ending index
% - num_sub: number of submatrices
% - n_x: system state
% - n_u: numeber of inputs
% - n_m: number of modes
% - inst: execution instant
% - labels: ideal labels of the past and future horizons points
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
    BPI_singles = NaN*ones(n_m+1,1);
    for i = 1 : n_m
        % N_h_i: number of elements associated with the mode_i in N_h
        N_h_i = sum(labels == i);
        % n_g: number of not null elements in g associated with the mode i
        n_g = 0;
        for j = 1 : num_sub
           if ind_H(j,3) == i 
              n_g = n_g + sum( abs(g(ind_H(j,1) : ind_H(j,2))) >= thr);
           end 
        end
        BPI_singles(i) = n_g/(n_u*N_h_i+n_x);
    end
    % Weights for multiple mode columns
    n_g = 0;
    % s: index for intervals
    s = 1;
    for j = 1 : num_sub
        % To avoid considering erroneously one column
        if s ~= ind_H(j,1)
           n_g = n_g + sum( abs(g(s : (ind_H(j,1)-1))) >= thr);
        end       
        % Updating index 
        s = ind_H(j,2) + 1; 
    end
    % Last columns (if is to avoid error)
    if s~= length(g)
        n_g = n_g + sum( abs(g(s : end)) >= thr);
    end
    BPI_singles(i+1) = n_g/n_x;
else 
    BPI_singles = NaN*ones(n_m+1,1);
end

end

