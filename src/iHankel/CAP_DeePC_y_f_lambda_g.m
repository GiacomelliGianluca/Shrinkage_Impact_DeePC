function [J_all, J_tr, J_g, g_opt, u_opt_seq, y_opt_seq] = CAP_DeePC_y_f_lambda_g(u_ini, y_ini, u_r, y_r, n_u, n_y, Q, R, N_h, N_ini, N_ex, H_u, H_y, D, lambda_g, n_m, ind_H, num_sub, slc_norm)

  
% FHCOP_rate
    % Solve the Finite Horizon Control Optimization Problem
    
    % Parameters
    % - u_ini: previous input
    % - y_ini: previous output
    % - u_r: The target input 
    % - y_r: The target output
    % - Q: The state weight
    % - R: The input weight
    % - N_h: The prediction horizon (number of steps)
    % - N_ini: Number of initialization steps
    % - H_u: Hankel matrix input
    % - H_y: Hankel matrix output
    % - curr: current iteration
    % - tot: total number of iterations
    % - n_m: number of modes
    % - num_sub: number of submatrices
    % - ind_H: indexes of the submatrices -Column 1: starting index -Column 2: ending index
    % - slc_norm: selected norm for CAP


    % Future Horizon
    L_h = N_h-N_ini;


    % Start cvx
    cvx_begin quiet
        % Declare the optimization variables
        variable u(n_u*L_h)
        variable y(n_y*L_h)
        variable g(D-N_h+1)

        % Cost function initialization
        J = 0;
        
        for ii=1:L_h

            % Cost function 
            J = J ...
                + (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)).' * Q * (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)) ...
                + (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1)).' * R * (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1));
        end

        CAP = 0;
        % Car: Cardinality of each group (included multiple mode columns)
        Car = zeros(n_m+1,1);
        % CAP_i: cap of mode i (included multiple mode columns)
        CAP_i = cell(n_m+1,1);
        for i = 1 : n_m
            % Initialization
            CAP_i{i,1} = [];
            for j = 1 : num_sub
               if ind_H(j,3) == i 
                  Car(i) = Car(i) + length(g(ind_H(j,1) : ind_H(j,2)));
                  CAP_i{i,1} = [CAP_i{i,1}; g(ind_H(j,1) : ind_H(j,2))];
               end 
            end
            CAP = CAP + sqrt(Car(i))*norm(CAP_i{i,1}, slc_norm);
        end
        % Weights for multiple mode columns
        % Initialization
        CAP_i{i+1,1} = [];
        % s: index for intervals
        s = 1;
        for j = 1 : num_sub
            % To avoid considering erroneously one column
            if s ~= ind_H(j,1)
               Car(i+1) = Car(i+1) + length(g(s : (ind_H(j,1)-1)));
               CAP_i{i+1,1} = [CAP_i{i+1,1} ; g(s : (ind_H(j,1)-1))];
            end       
            % Updating index 
            s = ind_H(j,2) + 1; 
        end
        % Last columns (if is to avoid error)
        if s~= length(g)
           Car(i+1) = Car(i+1) + length(g(s : end));
           CAP_i{i+1,1} = [CAP_i{i+1,1} ; g(s : end)];
        end

        CAP = CAP + sqrt(Car(i+1))*norm(CAP_i{i+1,1}, slc_norm);       


        % Adding the CAP at the cost function
        J = J + lambda_g*CAP;

        minimize(J)
        subject to
            % Model constraints
            H_u(1:N_ini,:)*g==u_ini;
            H_u(N_ini+1:end,:)*g==u;
            % H_y(1:N_ini,:)*g == y_ini + sigma_y;
            H_y(1:N_ini,:)*g == y_ini;
            H_y(N_ini+1:end,:)*g==y;
            % % Terminal constraints
            
            % Input constraints (same as hysel)
            -50 <= u <= 50


    cvx_end

    % Extract the cost function and the optimal weightings
    J_all = J;
    J_tr = 0;
    for ii=1:L_h
        % Cost function 
        J_tr = J_tr ...
             + (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)).' * Q * (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)) ...
             + (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1)).' * R * (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1));
    end
    J_g = 0;
    for i = 1 : (n_m+1)
        J_g = J_g + lambda_g*sum(sqrt(Car(i))*norm(CAP_i{i,1}, slc_norm));
    end
    % J_sig = norm(sigma_y,1);
    g_opt = g;
    % Extract the optimal input
    u_opt_seq = u;
    % Extract the optimal output prediction (for the sequence)
    y_opt_seq = y;



end