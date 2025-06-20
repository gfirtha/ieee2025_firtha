clear
%close all

rSSD = 2;         % SSD radius
nSSD = 360*2;       % secondary source number
xSource = [ 1, 0 ];    % source position
xRef = [0,0];
c = 343;
dx = 0.0025;
dx = 4*dx;

freq  = 1000;

phi = (0:nSSD-1)'/nSSD'*2*pi-pi;
dphi = mean(diff(phi));
x0 = rSSD*[cos(phi) sin(phi)];
n0 = [ -cos(phi) -sin(phi) ];


rhoP = sqrt(sum( (x0-xSource).^2,2 ));
rhoG = sqrt(sum(x0.^2,2));
kh = (x0-xSource)./rhoP;
khn = sum(n0 .* kh,2);
F = -(-1+2*all(khn<0)); % Focused flag
if F == 1
    win = double(khn>=0);
elseif F == -1
    ks = - (xRef-xSource)/norm(xRef-xSource);
    win = ((kh*ks')>0);
 %   win = (abs(kh*ks')).*((kh*ks')>0);
end
k0 = 2*pi*freq/c;

dref = rhoG;
Dx = win.*sqrt(8*pi*1i*k0).*sqrt(rhoP.*dref./(rhoP+F*dref)).*(khn).*exp(-F*1i*k0.*rhoP)./rhoP/4/pi;

%%

x = (-1.25*rSSD:dx:1.25*rSSD);
y = (-1.25*rSSD:dx:1.25*rSSD);
[X,Y] = meshgrid(x,y);
Xfield = [X(:),Y(:)];
field = zeros(size(Xfield,1),1);
R = sqrt(  (Xfield(:,1)-x0(:,1)').^2 + (Xfield(:,2)-x0(:,2)').^2  );

dl = mean(rhoG*2*pi/size(x0,1));
Psynth = reshape(  sum(Dx.'.*exp(-1i*k0*R)./R/4/pi*dl,2)  ,length(y),length(x));
R0 = sqrt( (X-xSource(1)).^2 + (Y-xSource(2)).^2  );
Pref = 1/(4*pi)*exp(-1i*k0*R0)./R0;
P0 = 1/(4*pi)*1./norm(xSource);

%%
% Vector from source to each grid point
D = Xfield - xSource;  % direction vectors
Dirx = D(:,1); Diry = D(:,2);
a = Dirx.^2 + Diry.^2;
b = 2*(xSource(1)*Dirx + xSource(2)*Diry);
c = xSource(1)^2 + xSource(2)^2 - rSSD^2;
disc = b.^2 - 4*a.*c;
valid = disc >= 0;
t1 = nan(size(disc));
t2 = nan(size(disc));
sqrt_disc = sqrt(disc(valid));
t1(valid) = (-b(valid) + sqrt_disc) ./ (2*a(valid));
t2(valid) = (-b(valid) - sqrt_disc) ./ (2*a(valid));
P1 = xSource + [Dirx.*t1, Diry.*t1];
P2 = xSource + [Dirx.*t2, Diry.*t2];
d1 = sqrt(sum((P1 - xSource).^2, 2));
d2 = sqrt(sum((P2 - xSource).^2, 2));
useP1 = d1 < d2;
x0Stat = P1;
x0Stat(~useP1,:) = P2(~useP1,:);

kvec = F*(x0Stat - xSource)./sqrt(sum((x0Stat - xSource).^2,2));
w0 = zeros(size(x0Stat, 1), 1);
for i = 1:size(x0Stat, 1)
    diffs = x0 - x0Stat(i, :);                % Mx2
    distsSquared = sum(diffs.^2, 2);          % Mx1
    [~, minIdx] = min(distsSquared);          % index of closest point
    w0(i) = win(minIdx);              % corresponding function value
    d0(i) = dref(minIdx);
end

xRefStat  = x0Stat + kvec.*d0';


A0 = w0.*sqrt( dist(xRefStat-x0Stat)./dist(xRefStat-xSource).*dist(Xfield-xSource)./dist(Xfield-x0Stat) );
A0 = reshape(A0, size(X));
ixs0 = X>xSource(1);
A0(ixs0)=A0(ixs0).*(-1i); 


R0 = sqrt( (X-xSource(1)).^2 + (Y-xSource(2)).^2  );
Pref = 1/(4*pi)*exp(-1i*k0*R0)./R0;
Pref(ixs0) = 1/(4*pi)*exp(1i*k0*R0(ixs0))./R0(ixs0);

%%
mask = sqrt( X.^2+Y.^2 )<=rSSD;
mask1 = double(mask);
mask1(mask==0) = nan;
mask2 = double(~mask);
mask2(mask2==0) = nan;

q= 10;
ftsize = 13;
scale_window = 0.25; % Adjust this value to scale window visualization size

pos = [ 0.18 0.18 .8 .8];
fig = figure('Units','points','Position',[150,150,200,200]);
p = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(Psynth).*mask1,'FaceAlpha',1);
hold on
pc2 = pcolor(x,y,real(Pref.*mask2),'FaceAlpha',0.5);
shading interp
axis equal tight
win_x = x0(:,1) + win .* n0(:,1) * scale_window;
win_y = x0(:,2) + win .* n0(:,2) * scale_window;
fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [245,154,159]/255, ...
     'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot(win_x, win_y, 'Color', [162,50,52]/255, 'LineWidth', 1.5);
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
clim([-1,1]*1/4/pi/norm(xSource));


% Arrow pointing from xSource to xRef
quiver(xSource(1), xSource(2), xRef(1)-xSource(1), xRef(2)-xSource(2), 0, ...
    'Color', 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off');

% Text annotation positioned slightly offset from the arrow midpoint
midpoint = (xSource + xRef) / 2;
text(midpoint(1)-0.8, midpoint(2), '$\mathbf{k}_{\mathrm{s}}$', ...
     'Interpreter', 'latex', 'FontSize', ftsize+2, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'k');


draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.05 );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'Psynth_circular_full' ,'-dpng')

%%

pos = [ 0.18 0.18 .8 .8];
fig = figure('Units','points','Position',[150,150,200,200]);
p = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(A0.*Pref).*mask1,'FaceAlpha',1);
hold on
pc2 = pcolor(x,y,real(Pref.*mask2),'FaceAlpha',0.5);
shading interp
axis equal tight
win_x = x0(:,1) + win .* n0(:,1) * scale_window;
win_y = x0(:,2) + win .* n0(:,2) * scale_window;
fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [245,154,159]/255, ...
     'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot(win_x, win_y, 'Color', [162,50,52]/255, 'LineWidth', 1.5);
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
clim([-1,1]*1/4/pi/norm(xSource));
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.05 );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'Psynth_circular_ideal' ,'-dpng')
%%

pos = [ 0.18 0.18 .8 .8];
fig = figure('Units','points','Position',[150,150,200,200]);
p = axes('Units','normalized','Position',pos(1,:));

pcolor(x,y,real(A0.*Pref-Psynth).*mask1,'FaceAlpha',1);
hold on
pc2 = pcolor(x,y,real(Pref.*mask2),'FaceAlpha',0.5);
shading interp
shading interp
axis equal tight
win_x = x0(:,1) + win .* n0(:,1) * scale_window;
win_y = x0(:,2) + win .* n0(:,2) * scale_window;
fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [245,154,159]/255, ...
     'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot(win_x, win_y, 'Color', [162,50,52]/255, 'LineWidth', 1.5);
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
clim([-1,1]*1/4/pi/norm(xSource));
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.05 );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'Psynth_circular_diffr' ,'-dpng')

function r = dist(x)

r = sqrt( sum( x.^2 , 2)  );
end

