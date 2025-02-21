function u_bar = computing_u_equil_PWA(y_bar, th)
% The function computes the corresponding input of a desired output for a PWA system

% y_bar: equilibrium outputs
% th: system parametes

% u_bar: equilibrium inputs

% Parameters
a1       =       th(1,1);     
a2       =       th(2,1);      
b1       =       th(3,1);      
b2       =       th(4,1);     

% Allocating u_bar
u_bar = NaN*ones(length(y_bar),1);

% Computing the equilibria
for i = 1 : length(y_bar)
    % Model equations
    if y_bar(i) >= 0
        % Region 2
       u_bar(i) = (1-a2)/b2 * y_bar(i);
    else
        % Region 1
       u_bar(i) = (1-a1)/b1 * y_bar(i);
    end
end 




end

