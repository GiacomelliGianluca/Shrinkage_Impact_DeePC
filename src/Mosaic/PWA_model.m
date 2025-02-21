function z_step = PWA_model(z, u, e, th)
% Function that compute the 1-step ahead simulation of the PWA model

% Parameters
a1       =       th(1,1);     
a2       =       th(2,1);      
b1       =       th(3,1);      
b2       =       th(4,1);     


% States
y               =       z(1,1); 

% Model equations
if y >= 0
    % Region 2
   y_step = a2*y + b2*u;
else
    % Region 1
   y_step = a1*y + b1*u; 
end

% Assigning the output
z_step = y_step;
end