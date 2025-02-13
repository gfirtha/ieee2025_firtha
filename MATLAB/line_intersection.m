function xIntersect = line_intersection(xs, xr, x0, x1)
% xs, xr, x0, x1 are each 2×1 vectors (or 1×2, but treat consistently)

% Direction vectors
A = xr - xs;   % Direction of line 1
B = x1 - x0;   % Direction of line 2

% Right-hand side
RHS = x0 - xs;

% Solve the linear system [A, -B] * [t; s] = RHS
% Build the 2x2 matrix
M = [ A(1), -B(1);
    A(2), -B(2) ];

% Solve for [t; s]
ts = M \ RHS;  % Equivalent to inv(M)*RHS but more efficient

t = ts(1);
% s = ts(2);  % If needed, can also retrieve s

% Intersection point
xIntersect = xs + t * A;
end