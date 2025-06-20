function [Phi_out,A_out] = getDerivatives(x0_in, xr_in, xs_in, F_in)

syms x_0 y_0 z_0 x_s y_s x_r y_r F real

% Definitions of r_S and r_G
r_S = sqrt((x_0 - x_s)^2 + (y_0 - y_s)^2 + z_0^2);
r_G = sqrt((x_r - x_0)^2 + (y_r - y_0)^2 + z_0^2);

% Definitions of Phi and A
Phi = -(F*r_S + r_G);
A = 1/(r_S^2 * r_G);

% Analytical derivatives w.r.t. z_0
dPhi = diff(Phi, z_0);
dA   = diff(A, z_0);
ddPhi = diff(Phi, z_0,2);
ddA   = diff(A, z_0,2);
dddPhi = diff(Phi, z_0,3);
dddA   = diff(A, z_0,3);
ddddPhi = diff(Phi, z_0,4);
ddddA   = diff(A, z_0,4);

Phi_fun      = matlabFunction(Phi, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r,F]);
A_fun        = matlabFunction(A, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r]);
dPhi_fun = matlabFunction(dPhi, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r,F]);
dA_fun   = matlabFunction(dA, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r]);
ddPhi_fun = matlabFunction(ddPhi, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r,F]);
ddA_fun   = matlabFunction(ddA, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r]);
dddPhi_fun = matlabFunction(dddPhi, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r,F]);
dddA_fun   = matlabFunction(dddA, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r]);
ddddPhi_fun = matlabFunction(ddddPhi, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r,F]);
ddddA_fun   = matlabFunction(ddddA, 'Vars',[x_0,y_0,z_0,x_s,y_s,x_r,y_r]);

params = [[x0_in(1:2),0] xs_in(1:2) xr_in(1:2) F_in ];
Phi_val      = Phi_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
A_val        = A_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7));
dPhi_val = dPhi_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
dA_val   = dA_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7));
ddPhi_val = ddPhi_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
ddA_val   = ddA_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7));
dddPhi_val = dddPhi_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
dddA_val   = dddA_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7));
ddddPhi_val = ddddPhi_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7),params(8));
ddddA_val   = ddddA_fun(params(1),params(2),params(3),params(4),params(5),params(6),params(7));


Phi_out = [Phi_val, dPhi_val, ddPhi_val, dddPhi_val, ddddPhi_val];
A_out = [A_val, dA_val, ddA_val, dddA_val, ddddA_val];
end

