function [xn, d, z, y] = temp_PWA_sys_sim(x, u, params)
% [xn, d, z, y] = temp_PWA_sys_sim(x, u, params)
% simulates the hybrid system one step ahead.
% Parameters:
%   x: current state
%   u: input
%   params: structure containing values for
%           all symbolic parameters
% Output:
%   xn: state in the next timestep
%   u: output
%   d, z: Boolean and real auxiliary variables
%
% HYSDEL 2.0.5 (Build: 20090715)
% Copyright (C) 1999-2002  Fabio D. Torrisi
% 
% HYSDEL comes with ABSOLUTELY NO WARRANTY;
% HYSDEL is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
if ~exist('x', 'var')
	error('error:  current state x not supplied');
end
x=x(:);
if ~all (size(x)==[1 1])
	error('error: state vector has wrong dimension');
end
if ~exist('u', 'var')
	error('error: input u not supplied');
end
u=u(:);
if ~all (size(u)==[1 1])
	error('error: input vector has wrong dimension');
end

d = zeros(1, 1);
z = zeros(1, 1);
xn = zeros(1, 1);
y = zeros(1, 1);

if (u(1) < -50) | (u(1) > 50)
	error('variable u is out of bounds');
end
if (x(1) < -15) | (x(1) > 15)
	error('variable x is out of bounds');
end

% sign_p = x >= 0;
within((0) - (x(1)), -15, 15, 24);
if (0) - (x(1)) <= 0
	d(1) = 1;
else
	d(1) = 0;
end

% z1 = {IF sign_p THEN a2 * x + b2 * u ELSE a1 * x + b1 * u};
if d(1)
	within(((0.9) * (x(1))) + ((0.15) * (u(1))), -21, 21, 26);
	z(1) = ((0.9) * (x(1))) + ((0.15) * (u(1)));
else
	within(((-0.3) * (x(1))) + ((1.4) * (u(1))), -74.5, 74.5, 26);
	z(1) = ((-0.3) * (x(1))) + ((1.4) * (u(1)));
end

% x = z1;
xn(1) = z(1);

% y = x;
y(1) = x(1);

xn=xn(:);
y=y(:);
z=z(:);
d=d(:);


function within(x, lo, hi, line)
 if x<lo | x>hi 
 error(['bounds violated at line ', num2str(line), ' in the hysdel source']); 
 end
