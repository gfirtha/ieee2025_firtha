function [x0Stat,x0StatX,x0StatY,xRefStat] = get_stat_pos(X,Xfield, xSource, rSSD)

% %%
% Vector from source to each grid point
D = Xfield - xSource;  % direction vectors
Dx = D(:,1); Dy = D(:,2);

% Quadratic coefficients for line-circle intersection
a = Dx.^2 + Dy.^2;
b = 2*(xSource(1)*Dx + xSource(2)*Dy);
c = xSource(1)^2 + xSource(2)^2 - rSSD^2;

% Discriminant
disc = b.^2 - 4*a.*c;

% Only keep valid intersections
valid = disc >= 0;
t1 = nan(size(disc));
t2 = nan(size(disc));

% Compute roots only where valid
sqrt_disc = sqrt(disc(valid));
t1(valid) = (-b(valid) + sqrt_disc) ./ (2*a(valid));
t2(valid) = (-b(valid) - sqrt_disc) ./ (2*a(valid));

% Evaluate points on line
P1 = xSource + [Dx.*t1, Dy.*t1];
P2 = xSource + [Dx.*t2, Dy.*t2];

% Choose the point closest to xSource
d1 = sqrt(sum((P1 - xSource).^2, 2));
d2 = sqrt(sum((P2 - xSource).^2, 2));
useP1 = d1 < d2;

% Final intersection point closest to xSource
x0Stat = P1;
x0Stat(~useP1,:) = P2(~useP1,:);

% Optionally reshape to grid
x0StatX = reshape(x0Stat(:,1), size(X));
x0StatY = reshape(x0Stat(:,2), size(X));

k0 = (x0Stat - xSource)./sqrt(sum((x0Stat - xSource).^2,2));
xRefStat  = x0Stat + k0*rSSD;



end

