function [H, index_H, f] = building_Hankel_w_index(Data, index, T_h, d, n_x, flag_PE)
%Function that constructs Hankel matrix

% Data: data retrieved from another vector
% index: index of the data in their original vector
% T_h: prediction horizon. It is composed by the past horizon T_ini and the
% future horizon L (T_h = T_ini + L). The initial instant assumed for the
% Hankel matrix is t = 0
% d: data dimensionality (corresponing to the number of elements for an
% istant t
% n_x: state dimensionality
% flag_PE: 1: test the persistency of excitation 0: no test


% H: Hankel matrix
% index_H: Hankel matrix realized with the indexes of the data

N_Data = length(Data);


%Partitioning raws and columns
% raw
raw_data = Data(1:d*T_h);
column_data = Data(d*T_h:N_Data);

% Same for the index
raw_index = index(1:d*T_h);
column_index = index(d*T_h:N_Data);

H = hankel(raw_data, column_data);
index_H = hankel(raw_index, column_index);

if flag_PE && (PE_check(Data, (d*(T_h + n_x))) ~= 1)
   keyboard
   H = [];
   index_H = [];
   % Subhankel not built
   f = 1;
else 
   % Subhankel built
   f = 0;
end



end