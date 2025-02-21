function obs_flag = obs_check(PWA_sys, obs_ind)
% obs_check checks the observability of the PWA system

% - PWA_sys: (raw) array of structs, where each struct is a system
%            having fields a, c
% - obs_ind: joint observability index


% n_sys: number of system
n_sys = length(PWA_sys(:,1));

% Extracting the elements 
A_PWA = cell(n_sys,1);
C_PWA = cell(n_sys,1);

for i = 1 :n_sys
    A_PWA{i,1} = PWA_sys(i).a;
    C_PWA{i,1} = PWA_sys(i).c;
end

% Extracting the system order 
n_x = length(A_PWA{i,1});

% Building the joint observability matrix
Obs_1 = [];
Obs_2 = [];
for i = 1 : obs_ind
    if i == 1
        Obs_1 = [Obs_1; C_PWA{1,1}];
        Obs_2 = [Obs_2; C_PWA{2,1}];
    else
        Obs_1 = [Obs_1; Obs_1(i-1,:)*A_PWA{1,1}];
        Obs_2 = [Obs_2; Obs_2(i-1,:)*A_PWA{2,1}];
    end     
end
% Obs_12 joint observability matrix
Obs_12 = [Obs_1 Obs_2];

% no_obs:  flag regarding observability (if no_obs = 1, the PWA system is
% not observable)
no_obs = 0;

% Checking condition 1
if rank(Obs_12) ~= 2*n_x
    % The PWA system is not observable
    keyboard;
    no_obs = 1;
end

% Checking condition 2
if  (rank((Obs_1 - Obs_2)*A_PWA{1,1}) ~= n_x) || (rank((Obs_2 - Obs_1)*A_PWA{2,1}) ~= n_x)
    % The PWA system is not observable
    keyboard;
    no_obs = 1;
end

% Checking condition 3 for mode 1
if  (rank(A_PWA{1,1} - A_PWA{2,1}) ~= n_x) || (rank(A_PWA{1,1}^2 - A_PWA{2,1}^2) ~= n_x) || (rank(A_PWA{1,1}^2 - A_PWA{1,1} * A_PWA{2,1}) ~= n_x)
    % The PWA system is not observable
    keyboard;
    no_obs = 1;
end

% Checking condition 3 for mode 2
if  (rank(A_PWA{2,1} - A_PWA{1,1}) ~= n_x) || (rank(A_PWA{2,1}^2 - A_PWA{1,1}^2) ~= n_x) || (rank(A_PWA{2,1}^2 - A_PWA{2,1} * A_PWA{1,1}) ~= n_x)
    % The PWA system is not observable
    keyboard;
    no_obs = 1;
end

if no_obs == 1
    obs_flag = 0;
else 
    obs_flag = 1;
end

end

