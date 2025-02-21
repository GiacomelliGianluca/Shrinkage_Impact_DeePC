% Main script for an interpretable DeePC for PWA systems using a Mosaic data-driven model

close all 
clear 
clc

% plt_set: Plot settings
plt_set.y_font_dim = 16;
plt_set.x_font_dim = 11;
plt_set.title_font_dim = 14;
plt_set.thick = 1.6;
plt_set.plot_dim_x = 8; 
plt_set.plot_dim_y = 6; 
plt_set.plot_unit = 'centimeters'; 
plt_set.fontname = 'Times';
plt_set.fontsize = 12;

% Dataset directory
dataset_dir = 'Set data path';

%% System 

% n_u Input Size
n_u = 1;
% n_y Output Size
n_y = 1;
% n_x State Size
n_x = n_y;

% Definind the system as HYSDEL model
S = mld('PWA_sys');
% Retrieving the associated PWA system
P = pwa(S);
% System parameters
sys1.a = -0.3;
sys2.a = 0.9;
sys1.b = 1.4;
sys2.b = 0.15;
sys1.c = 1;
sys2.c = 1;
% Parameter vector
th       =       [sys1.a;sys2.a;sys1.b;sys2.b];

% joint_index: joint observability index
joint_index = joint_obs_index([sys1; sys2]);
% Observability check (obs_flag = 1: yes; obs_flag = 0: no)
obs_flag = obs_check([sys1; sys2], joint_index);

%% Reading the data

% Selecting the datasest to build the Hankel
string_data_hankel = 'Data_two_modes.mat';
Data_u_model = cell2mat(struct2cell(load([dataset_dir,string_data_hankel], "Data_u")));
Data_y_model = cell2mat(struct2cell(load([dataset_dir,string_data_hankel], "Data_y")));
step_model = cell2mat(struct2cell(load([dataset_dir,string_data_hankel], "step")));

% Dataset length (Assuming that Data_u and Data_y have the same data)
N_D = length(Data_u_model);

%% Parameters of DDPC 

% Number of initialization steps (past horizon)
N_ini = 25;
% Future horizon (where L = 1 corresponds with t=0)
L  = 19;
% Prediction Horizon
N_h = N_ini + L;
% Shrinkage regularization 
lambda_g = 10;
% Control Horizon
T_c = 50;
% slc_norm: selected norm for the penalty (only iCAP-DeePC)
slc_norm = 2;

%% Partition the data

% Forming the data vector
Data = [Data_u_model' Data_y_model'];

% Initial data plot
Intial_data = figure(1);
plot(Data(:,1),Data(:,2),'k*','MarkerSize',5);
grid on,
str_x = '$u_{d_k}$';
xlabel(str_x,'Interpreter','latex');
str_y = '$y_{d_k}$';
ylabel(str_y,'Interpreter','latex');
str_title = 'Dataset';
title(str_title, 'Interpreter','latex', 'FontWeight', 'bold');

% Forming the data vector for the N_ini clustering
Data_cl = NaN*ones(N_D-N_ini, 2*N_ini+1);
j = 1;
for i = 1 : N_D-N_ini
   Data_cl(j,:) = [Data_u_model(i:i+N_ini) Data_y_model(i:i+N_ini-1)];
   j = j + 1;
end


% n_m: Number of modes (it corresponds with the number of researched
% clusters)
n_m = 2;

% Clustering - ideal regions
[label_id, Centr_id] = Cluster_manual(Data_y_model, 0);

% kmeans - N_ini 
% [label_N_ini, Centr_N_ini] = kmeans(Data_cl,n_m, 'EmptyAction','error', OnlinePhase='on');

% Centroids that provide the results reported in the article
Centr_square_random_N_ini = [2.77905473965287	2.63928793947486	2.55697819314642	2.51571962616823	2.16404361370717	1.74057943925234	1.25378282153983	0.684712950600801	0.0504165554072096	-0.649106364040944	-1.41385580774366	-2.24062750333778	-3.07526924788607	-3.90182376502002	-4.17200445037827	-4.37926924788606	-4.52361815754339	-4.59659546061415	-4.60702892745883	-4.55491855807743	-4.44026435246996	-4.26306631063640	-3.79813529149978	-3.27248509123276	-2.82418869603917	-2.42801246105919	3.69218691588785	3.86093457943925	3.95028037383178	3.96022429906542	3.89076635514019	3.74190654205608	3.51364485981309	3.20598130841122	2.81607476635514	2.34392523364486	1.78953271028038	1.15289719626168	0.438504672897197	-0.279177570093458	-1.00014953271028	-1.64216822429907	-2.20523364485981	-2.68934579439252	-3.09450467289720	-3.42071028037383	-3.66796261682243	-3.83626168224299	-3.92560747663552	-3.93600000000000	-3.86743925233645
                            -4.83377597402597	-4.62660714285714	-4.48977272727273	-4.40381168831169	-3.94137012987013	-3.39259956709957	-2.76778138528139	-2.06691558441558	-1.29000216450216	-0.437041125541125	0.491967532467533	1.49312770562771	2.48706709956710	3.46850649350649	3.77077922077922	3.99388528138528	4.13782467532468	4.20259740259740	4.18820346320346	4.09464285714286	3.92191558441559	3.67002164502165	3.09545454545455	2.44696969696970	1.89245129870130	1.40121753246753	-4.59672727272727	-4.77227272727273	-4.84954545454545	-4.83081818181818	-4.71609090909091	-4.50536363636364	-4.19863636363636	-3.79590909090909	-3.29718181818182	-2.70245454545454	-2.01172727272727	-1.22500000000000	-0.347727272727273	0.529545454545454	1.40681818181818	2.18409090909091	2.86136363636364	3.43863636363637	3.91590909090909	4.29318181818182	4.57045454545455	4.74772727272727	4.82500000000000	4.80681818181819	4.69318181818182];
[label_N_ini, Centr_N_ini,~, dist_point_N_ini] = kmeans(Data_cl,n_m, 'EmptyAction','error','start', Centr_square_random_N_ini, OnlinePhase='on');
% Adding the first unclassified rho traj samples (labeled with -1)
label_N_ini = [-1*ones(N_ini,1); label_N_ini];


% Cluster data plot - knowing regions
data1_string = '$u_{d_k}$';
data2_string = '$y_{d_k}$';
title_string = 'Partitioned dataset - Ideal';
Clustered_data_id = Cluster_data_plot_manual(Data, label_id, title_string, data1_string, data2_string, 2);

% Cluster data plot - k-means  
data1_string = '$u_{d_k}$';
data2_string = '$y_{d_k}$';
title_string = ['Partitioned dataset - ', num2str(2*N_ini),' features per data point'];
Clustered_data = Cluster_data_plot_manual(Data, label_N_ini, title_string, data1_string, data2_string, 3);

% f_inv: flag to keep track that the labels are inverted (relevant with
% k-means clustering)
f_inv = input("Are the labels inverted? 1: Yes, 0: No ");

% State the cluster type (id: ideal, cl: kmeans)
% str_cl_type = 'id';
str_cl_type = 'cl';
% Selecting the scenario
% % Ideal clustering (knowing the regions)
% Centr = Centr_id;
% label = label_id;
% k-means N_ini
Centr = Centr_N_ini;
label = label_N_ini;
% Updating the dataset (unlabeled data are discarded)
N_D = N_D - N_ini;

%% DDPC 

% Number of simulation steps
N = N_ini + T_c;
% Execution period steps (1: we execute at each step)
N_ex = 1;

% Setting the references 
% Output reference evolving from x_0 to x_final based on spline 
% Initial point of the trajectory
y_0_r_ini = -10;
% Final point of the trajectory
y_f_r_ini = -10;
% Defining the middle points for spline (spline since it is cubic needs 4 points to work)
if y_0_r_ini < y_f_r_ini
    y_mid_1_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.25 + y_0_r_ini;
    y_mid_2_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.75 + y_0_r_ini;
else
    y_mid_1_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.75 + y_f_r_ini;
    y_mid_2_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.25 + y_f_r_ini;
end
% Defining the points for spline
y_spline = [ y_0_r_ini  y_mid_1_r_ini  y_mid_2_r_ini   y_f_r_ini];
x_spline = [1 N_ini/4  3/4*N_ini N_ini];
% Defining the reference vector
y_r_ini = spline(x_spline, y_spline, 1:1:N_ini);


% Periods of the reference signal
T_r_1 = round(L/2)-1; % Before the switch 
T_r_2 = T_c - T_r_1; % After the switch 

% Future horizon references 
% Output reference sequence  
y_r_1 = -10*ones(1,T_r_1);
y_r_2 = 10*ones(1, T_r_2);
% y_r_ch: output reference control horizon
y_r_ch = [y_r_1     y_r_2];
% Checking the correct length
if length(y_r_ch)~=T_c
    keyboard;
end
y_r_f = [y_r_ch  y_r_ch(end)*ones(1,L)]';

% Input reference 
u_r_f = computing_u_equil_PWA(y_r_f, th);

str_xaxis = '$t$ [samples]';
% Steps (x-axis)
steps = 0: (T_c - 1);

%% Defining the MPC for the initial trajectory

% Signals in cost function
refs.y = 1;   % output references (element index 1)
% refs.u = 1;   % input references (element index 1)
% Weighting matrices cost function
Q.y = 1e3*eye(n_y);      % weight cost function output
% Q.u = 1*eye(n_u);      % weight cost functioon input
Q.rho = Inf;           % hard constraints (no slack)
% Type of norm of cost function elements
Q.norm = 2;         

% Constraints (Same as the ones in the hysel model)
limits.umin = -50;
limits.umax = 50;
limits.xmin = -15;
limits.xmax = 15;
% % Solver
mipsolver = 'gurobi';

% C: MPC Controller
C = hybcon(S,Q,N_h,limits,refs,mipsolver);

% Weights for the DeePC 
% Q Matrix DeePC (weight of y)
Q = 1*eye(n_y);
% R Matrix DeePC (weight of u)
R = 1*eye(n_u);

% H_u: (mosaic) Input Hankel Matrix
% index_M: matrix of the indexes of the data of the mosaic Hankel matrix
% centr_M: matrix of the centroids (modes) associated with the data of the mosaic Hankel matrix
% M_label: matrix of the same size of H_M, containing the labels of the partitioned data
[H_u, index_M, M_label, ~, n_m_u] = building_Mosaic_PWA(Data_u_model', N_h, n_m, label, Centr, n_u, n_x, 1);
% H_u: (mosaic) Output Hankel Marix
[H_y, ~, ~, centr_M, n_m_y, label_del] = building_Mosaic_PWA(Data_y_model', N_h, n_m, label, Centr, n_y, n_x, 0);

% Checking if the number of availiable matrices is the same
if n_m_u ~= n_m_y
    keyboard;
else
    n_m = n_m_y;
end
% Checking if some data has been deleted
if ~isempty(label_del)
    % Updating the Centroids
    Centr(label_del,:) = [];
    % Identify the delated data 
    Data_del = Data_y_model(label == label_del);
    % Updating the number of data
    N_D = N_D - length(Data_del);
end

% Declaring the switching condition
cond.y = 0;

% Plotting the Mosaic matrix
next_fig_index = 4;
fig_map = Hankel_map_PWA_m(str_cl_type, n_m, M_label, plt_set, next_fig_index, 'Mosaic');

% Plotting the data strcuture higlighting the missclassified points
% M_label_miss: matrix containing the labels of data and the
% misclassfication indicated as 0
[M_label_miss, fig_map_miss] = Hankel_map_miss_PWA(str_cl_type, fig_map, label, label_id, f_inv, M_label, index_M, N_ini, plt_set, 5, 'Mosaic');


% Plotting the number and avarage magnitude (H_miss_mag_avg) of missclassified data for a column
% M_miss_num: number of misclassification 
% M_miss_mag_avg_u: input average magnitude for each column 
% M_miss_mag_avg_y: output average magnitude for each column
[M_miss_num, M_miss_mag_avg_u, M_miss_mag_avg_y, fig_col_miss_num, fig_col_miss_mag_avg] = DS_miss_num_mag_PWA(str_cl_type, M_label_miss, H_u, H_y, plt_set, 6, 7);


% Identify the submatrixes for each mode
% sub_num: number of submatrices (it corresponds with the number of clusters)
sub_num = n_m;
% sub_ind_H: matrix containing the indexes in the hankel matrix of the submatrices 
% fig_H_sub: figure handle of the submatrices plot
[sub_ind_H, fig_H_sub] = Hankel_sub_Mos(M_label, N_h, n_x, n_u, n_m, fig_map, plt_set, 8);
% Identify submatrices for each mode without missclassifications 
% sub_ind_H_miss: matrix containing the indexes in the hankel matrix of the submatrices (if any)
% sub_num_miss: number of submatrices
% fig_H_sub_miss: figure handle of the submatrices plot
[~, sub_ind_H_miss, sub_num_miss, fig_H_sub_miss] = Hankel_sub_miss_PWA(str_cl_type, M_label_miss, sub_ind_H, sub_num, H_u, N_h, n_x, n_u, fig_map_miss, plt_set, 9);

% Allocation input applied (First N_ini steps is determined by an MPC, then
% the remaining L steps is decided through DDPC)
u = NaN*ones(n_u*(N_ini+T_c),1);
% Allocation system output (First N_ini steps is determined by an MPC, then
% the remaining L steps is decided through DDPC)
y = NaN*ones(n_y*(N_ini+T_c),1);
% Allocation overall cost
J_opt = NaN*ones(T_c,1);
% Allocation tracking cost
J_tr = NaN*ones(T_c,1);
% Allocation reg on g cost
J_g = NaN*ones(T_c,1);

% Allocation of the optimal inpute sequence (computed by the DDPC)
u_opt_seq = NaN*ones(n_u*L,T_c);
% Allocation of the predicted output (computed by the DDPC)
y_pre_seq = NaN*ones(n_y*L,T_c);
% Allocation of optimal weights 
g_opt = NaN*ones(N_D-n_m*(N_h-1),T_c);

%% Metrics computation 
% RMSE_u (Prediction horizon) 
RMSE_u = NaN;
% RMSE_y (Prediction horizon) 
RMSE_y = NaN;
% BPI_inst_single 
BPI_inst_single = cell(n_m,1);
% Initializing cell
for i = 1 :n_m
    % For every mode, in the cell is allocated a matrix corresponding with rows as the instants in L and the columns as hyperparameter
    % values.
    BPI_inst_single{i,1} = NaN*ones(T_c, 1);
end

% thr_null: value under which we assume null a component of g 
thr_null = 1e-2;


%% Execution

% Allocation of the vector of the simulation (each column corresponds to a state)
y_sim = zeros(N, n_y);
% Initialization IC ode_45
y_sim (1, 1) = y_r_ini(1,1);

% Reallocation of u_ini, y_ini, and labels for the past and future horizons for BPI calculation (-1: unassigned)
u_ini = zeros(n_u*N_ini, 1);
y_ini = NaN*ones(n_y*N_ini, 1);
% Labels_id (ideal, model-based) 
Labels_id = -1*ones(N_h,1);
% Labels_k (k-means) 
Labels_k = -1*ones(N_h,1);
for k= 1 : (N_ini + T_c)

    % Past horizon
    if k <= N_ini

        % MPC case 
        if k == 1
            % Setting the references
            r.y = y_r_ini';
            % MPC for the trajectory of the past horizon (control period N_ini)
            [~, u_c_seq, ~, ~, TT, y_c_seq] = sim(C, S, r, y_sim(k, :)', N_ini);
        end
         % Update MPC control actions
         u(k) = u_c_seq(((k-1)*n_u+1):n_u*k,1);
     end


     % Computation of the optimal control action
     if k >= N_ini + 1          % Checking if the DeePC has to be executed
         % Index control period
         fut = k - N_ini;
         % iLasso-DeePC 
         [J_opt(fut), J_tr(fut), J_g(fut), g_opt(:,fut), u_opt_seq(:, fut), y_pre_seq(:, fut)] = DeePC_Mos_y_f_lambda_g(u_ini, y_ini, u_r_f(fut: fut + L - 1), y_r_f(fut: fut + L - 1), n_u, n_y, Q, R, N_h, N_ini, N_ex, H_u, H_y, N_D, n_m, lambda_g);
         % % iCAP-DeePC 
         % [J_opt(fut), J_tr(fut), J_g(fut), g_opt(:,fut), u_opt_seq(:, fut), y_pre_seq(:, fut)] = CAP_DeePC_Mos_y_f_lambda_g(u_ini, y_ini, u_r_f(fut: fut + L - 1), y_r_f(fut: fut + L - 1), n_u, n_y, Q, R, N_h, N_ini, N_ex, H_u, H_y, N_D, n_m, lambda_g, sub_ind_H, sub_num, slc_norm);
         % Update control action
         u(k) = u_opt_seq(1, fut);

         % Computation of the metrics on the selector
        if~isnan(y_pre_seq(:, fut))
            % Number of missclassifcations by kmeans 
            [~, Labels_id, Labels_k] = Num_miss(str_cl_type, [ u_ini y_ini ; u_r_f(fut: fut + L - 1) y_r_f(fut: fut + L - 1)], cond, Centr, Labels_k, fut, N_ini, f_inv);
            if sub_num > 0
                % BPI    
                BPI_singles = BPI_comp_m('id', g_opt(:,fut)', sub_ind_H, n_x, n_u, n_m, Labels_id, fut, N_ini, f_inv);         
                for f = 1 : n_m
                    BPI_inst_single{f,1}(fut , 1) = BPI_singles(f);
                end
            end
        else
            keyboard;
        end
     end 
      % System evolution
      y_sim (k+1, :) = PWA_model(y_sim (k, :),  u(k) , 0, th);
        
        
      % Update of the initial traj
      u_ini = [u_ini(n_u*2:end,1);    u(k)];
      y_ini = [y_ini(n_u*2:end,1);    y_sim(k, 1)];
    
     
 end

 % Saving the simulation output
 y = y_sim;

 if~isnan(y_pre_seq)
    % Metrics computation 
    % RMSE_u input (Prediction horizon) 
    RMSE_u = sqrt(1/T_c*sum((u(N_ini+1 : end, 1) - u_r_f(1: T_c)).^2));        
    % RMSE_y output (Prediction horizon) 
    RMSE_y = sqrt(1/T_c*sum((y(N_ini+1 : end-1, 1) - y_r_f(1: T_c)).^2));
end


%% Plotting the results 

% Input - min tracking error
fig_input = figure(10);
plot(steps, u_r_f(1:T_c),'k', 'linewidth',1), hold on, grid on,
plot(steps, u(N_ini+1:N_ini+T_c),'Color', [0.33 0.15 0.56],'linewidth',plt_set.thick), grid on,
xlim([0 49]);
str = str_xaxis;
x_title = xlabel(str,'Interpreter','latex');
str_y = '$u_t$';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Input tracking','Interpreter','latex', 'FontWeight', 'bold');
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_input.Units               = plt_set.plot_unit;
fig_input.Position(3)         = plt_set.plot_dim_x;
fig_input.Position(4)         = plt_set.plot_dim_y;
set(fig_input.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_input.PaperPositionMode   = 'auto';

% Output - min tracking error
fig_output = figure(11);
plot(steps, y_r_f(1:T_c),'k','linewidth',1), hold on, grid on,
plot(steps, y(N_ini+1 : end-1, 1),'Color', [0.33 0.15 0.56],'linewidth',plt_set.thick), 
str = str_xaxis;
x_title = xlabel(str,'Interpreter','latex');
xlim([0 49]);
str_y = '$y_t$';
y_title = ylabel(str_y,'Interpreter','latex');
str_title = 'Output tracking';
h_title = title(str_title, 'Interpreter','latex', 'FontWeight', 'bold');set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_output.Units               = plt_set.plot_unit;
fig_output.Position(3)         = plt_set.plot_dim_x;
fig_output.Position(4)         = plt_set.plot_dim_y;
set(fig_output.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_output.PaperPositionMode   = 'auto';


% Plotting the weights of the trajectories
fig_weights_all = plotting_weights_L(g_opt, sub_num, sub_ind_H, plt_set, 12); 
% Plotting the BPI of each mode along the prediction horizon
fig_BPI_single = plotting_BPI_single_m(BPI_inst_single, plt_set, 13);


%% Saving performances in txt
writetable(table(RMSE_u, RMSE_y),'Results.txt');
open('Results.txt');
