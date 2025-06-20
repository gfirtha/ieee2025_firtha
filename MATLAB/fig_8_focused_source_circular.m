clear
%close all

rSSD = 2;         % SSD radius
nSSD = 360*2;       % secondary source number
xSource = [ 1, 0 ];    % source position
xRef = [0,0];
c = 343;
dx = 0.0025;
dx = 4*dx;

freq  = [400,1000,10000];

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
    win1 = double(khn>=0);
elseif F == -1
    ks = - (xRef-xSource)/norm(xRef-xSource);
    win1 = ((kh*ks')>0);
    %   win = (abs(kh*ks')).*((kh*ks')>0);
end
k0 = 2*pi*freq/c;
[~,max_ix] = max(abs(khn.*win1));
win0 = tukeywin(length(find(win1==1)),0.5);
win0 = [win0;zeros(length(find(win1~=1)),1)];
win0 = circshift(win0,max_ix-round(length(find(win1==1))/2)-1);
win = win1.*win0;

dref = rhoG;

x = (-1.25*rSSD:dx:1.25*rSSD);
y = (-1.25*rSSD:dx:1.25*rSSD);
[X,Y] = meshgrid(x,y);
Xfield = [X(:),Y(:)];
field = zeros(size(Xfield,1),1);
R = sqrt(  (Xfield(:,1)-x0(:,1)').^2 + (Xfield(:,2)-x0(:,2)').^2  );
dl = mean(rhoG*2*pi/size(x0,1));

for fi = 1 : length(freq)
    k0 = 2*pi*freq(fi)/c;
    Dx = win.*sqrt(8*pi*1i*k0).*sqrt(rhoP.*dref./(rhoP+F*dref)).*(khn).*exp(-F*1i*k0.*rhoP)./rhoP/4/pi;
    Psynth{fi} = reshape(  sum(Dx.'.*exp(-1i*k0*R)./R/4/pi*dl,2)  ,length(y),length(x));
    R0 = sqrt( (X-xSource(1)).^2 + (Y-xSource(2)).^2  );
    Pref{fi} = 1/(4*pi)*exp(-1i*k0*R0)./R0;
end
P0 = 1/(4*pi)*1./norm(xSource);
%%
% freq_vec = logspace(log10(20),log10(5e3),256);
% k_vec = 2*pi*freq_vec/c;
% ddPhi = -(F./rhoP + 1./dref);
% SPA_corr = sqrt(2*pi./(k_vec.*abs(ddPhi))).*exp(1i*pi/4*sign(ddPhi));
% Dx_mx = -2*1i.*win.*k_vec.*SPA_corr.*(-khn)./rhoP.*exp(-F*1i*k_vec.*rhoP)/4/pi;
% G0 = 1/4/pi*exp(-1i*k_vec.*vecnorm(x0-xRef,2,2))./vecnorm(x0-xRef,2,2);
% P_transfer = sum( Dx_mx.*G0*dl , 1).'./(1/4/pi*exp(-1i*k_vec.'.*norm(xRef-xSource))./norm(xRef-xSource));

%%
scale_window = 0.25; % Adjust this value to scale window visualization size

mask = sqrt( X.^2+Y.^2 )<=rSSD;
mask1 = double(mask);
mask1(mask==0) = nan;
mask2 = double(~mask);
mask2(mask2==0) = nan;

unitIxs1 = find(win1==1);
unitIxs2 = find(win0==1);
cricticalIxs = [ unitIxs1(1), unitIxs1(end), unitIxs2(1),unitIxs2(end)];

LW = 1;
pos = [ 0.15 0.175 0.8 0.8 ];
ftsize = 13;

fig = figure('Units','points','Position',[150,150,250,250]);
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(Psynth{2}.*mask1),'FaceAlpha',1);
hold on
pcolor(x,y,real(Pref{2}.*mask2),'FaceAlpha',0.25);
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

xl = gca().XLim;
yl = gca().YLim;
for ix = 1:length(cricticalIxs)
    i = cricticalIxs(ix);
    p0 = x0(i,:);
    v = -kh(i,:);  % Direction vector (inward from x0 toward center or source)
    a = dot(v, v);
    b = 2 * dot(p0, v);
    c = dot(p0, p0) - rSSD^2;
    t_intersect = (-b + sqrt(b^2 - 4*a*c)) / (2*a); % Positive root
    p_end = p0 + t_intersect * v;
    if ix <= 2
        lstyle = ':';
    else
        lstyle = '--';
    end
    line([p0(1), p_end(1)], [p0(2), p_end(2)], ...
        'Color', 'black', 'LineWidth', LW, 'LineStyle', lstyle);
end
xlim(xl);
ylim(yl);

text(-1,0,'𝔍','FontSize',14)  % Represents the fraktur J
text(1.3,0, '𝔖', 'FontSize', 14)
text(0.5, 1.3, '𝔗_1', 'FontSize', 14)
text(0.5, -1.3, '𝔗_2', 'FontSize', 14)

q= 10;
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.05 );

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', 'circular_WFS_setup_foc' ,'-dpng')

%%
for n = 1: length(Psynth)
    pos = [ 0.17 0.2 .8 .78];
    ftsize = 13;
    fig = figure('Units','points','Position',[150,150,200,200]);
    p3 = axes('Units','normalized','Position',pos(1,:));

    pcolor(x,y,abs(20*log10(abs(Psynth{n})./abs(Pref{n}))))
    shading interp
    axis equal tight
    clim([0,8])
    hold on
    xl = gca().XLim;
yl = gca().YLim;
for ix = 1:length(cricticalIxs)
    i = cricticalIxs(ix);
    p0 = x0(i,:);
    v = -kh(i,:);  % Direction vector (inward from x0 toward center or source)
    a = dot(v, v);
    b = 2 * dot(p0, v);
    c = dot(p0, p0) - rSSD^2;
    t_intersect = (-b + sqrt(b^2 - 4*a*c)) / (2*a); % Positive root
    p_end = p0 + t_intersect * v;
    if ix <= 2
        lstyle = ':';
    else
        lstyle = '--';
    end
    line([p0(1), p_end(1)], [p0(2), p_end(2)], ...
        'Color', 'black', 'LineWidth', LW, 'LineStyle', lstyle);
end
xlim(xl);
ylim(yl);
    win_x = x0(:,1) + win .* n0(:,1) * scale_window;
    win_y = x0(:,2) + win .* n0(:,2) * scale_window;
    fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [245,154,159]/255, ...
        'FaceAlpha', 0.8, 'EdgeColor', 'none');
    plot(win_x, win_y, 'Color', [162,50,52]/255, 'LineWidth', 1.5);

    xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
    ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
    q= 10;
    draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.05 );
    set(gca,'FontName','Times New Roman');
    allAxesInFigure = findall(fig,'type','axes');
    set(allAxesInFigure,'FontSize',ftsize);
    set(gcf,'PaperPositionMode','auto');
    print( '-r300', sprintf('circular_WFS_foc_dA_mes_%i',n) ,'-dpng')
end

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
%%
pos = [ 0.15 0.175 0.8 0.8 ];
ftsize = 13;

fig = figure('Units','points','Position',[150,150,300,250]);
p2 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,abs(20*log10(abs(A0))))
shading interp
axis equal tight
clim([0,8])
c = colorbar;
c.Label.String = '[dB]';
hold on
win_x = x0(:,1) + win .* n0(:,1) * scale_window;
win_y = x0(:,2) + win .* n0(:,2) * scale_window;
fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [245,154,159]/255, ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot(win_x, win_y, 'Color', [162,50,52]/255, 'LineWidth', 1.5);


xl = gca().XLim;
yl = gca().YLim;
for ix = 1:length(cricticalIxs)
    i = cricticalIxs(ix);
    p0 = x0(i,:);
    v = -kh(i,:);  % Direction vector (inward from x0 toward center or source)
    a = dot(v, v);
    b = 2 * dot(p0, v);
    c = dot(p0, p0) - rSSD^2;
    t_intersect = (-b + sqrt(b^2 - 4*a*c)) / (2*a); % Positive root
    p_end = p0 + t_intersect * v;
    if ix <= 2
        lstyle = ':';
    else
        lstyle = '--';
    end
    line([p0(1), p_end(1)], [p0(2), p_end(2)], ...
        'Color', 'black', 'LineWidth', LW, 'LineStyle', lstyle);
end
xlim(xl);
ylim(yl);

q= 10;
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.05 );
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', 'circular_WFS_dA_est_foc' ,'-dpng')

function r = dist(x)

r = sqrt( sum( x.^2 , 2)  );
end

