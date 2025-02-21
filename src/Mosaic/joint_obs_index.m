function index_s = joint_obs_index(PWA_sys)
% joint_obs_index computes the joint observability index for a PWA system

% - PWA_sys: (raw) array of structs, where each struct is a system
%            having fields a, c

% n_sys: number of system
n_sys = length(PWA_sys(:,1));

% Extracting the elements 
A_PWA = cell(n_sys,1);
C_PWA = cell(n_sys,1);

for i = 1 :n_sys
    A_PWA{i,1} = PWA_sys(i).a;
    C_PWA{i,1} = PWA_sys(i).c;
end

% Computing the joint observability index
% Maximum index to inspect
index_max = 10;

index_s = 0;

Obs_i = [];
Obs_j = [];

for i = 1 : index_max
    if i == 1
        Obs_i = [Obs_i; C_PWA{1,1}];
        Obs_j = [Obs_j; C_PWA{2,1}];
    else
        Obs_i = [Obs_i; Obs_i(i-1,:)*A_PWA{1,1}];
        Obs_j = [Obs_j; Obs_j(i-1,:)*A_PWA{2,1}];
    end
    index = rank([Obs_i Obs_j]);

    if index > index_s
        index_s = index;
        % c: counter of the number of times the rank is not increasing
        c = 0;
    else 
        c = c + 1;
    end    
end

end

