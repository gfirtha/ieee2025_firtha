clear
close all
addpath(genpath( 'Functions\' ))

Lssd = 1.5;                 % Half-length of integration domain in x-direction
dxSSD = 1e-3;
dxMesh = 1e-2;
dxDif = 0.001;
taper = 0.5;
freq = logspace(-1,log10(10e3),256)'; % Frequency vector (logarithmic scale)
c = 343;                    % Speed of sound in air (m/s)
f1  = 400;
f2  = 1000;
f3  = 10000;
xRef = [0,1];       % Receiver/reference position
xSource = [0,-1];        % Source position
phiLim = pi/4;

x0 = (-Lssd:dxSSD:Lssd)';        % x-axis spatial vector
y0 = 0;                   % Fixed y-position of the integration line
win = tukeywin(length(x0),taper); % Tukey window applied to x-axis

%% Calculate radiated field
F = sign(y0-xSource(2));
aEvanComp = (xRef(2)-F*xSource(2)) / (xRef(2)-xSource(2));
rS = sqrt( (x0-xSource(1)).^2 + (y0-xSource(2)).^2 );

dref = 1./(abs(xRef(2)./(xRef(2)-xSource(2))))./rS;
SPA_corr = sqrt(2*pi./abs(dref)).*exp(-1i*pi/4*sign(dref));
Dx0_prop = 2/(4*pi)*win.*SPA_corr.*(y0-xSource(2))./rS.^2;

k1 = 2*pi*f1 / c;
k2 = 2*pi*f2 / c;
k3 = 2*pi*f3 / c;
Lx = 2;
Ly = 2;
x = (-Lx:dxMesh:Lx)';
y = (dxMesh:dxMesh:Ly)';
[X,Y] = meshgrid(x,y);
q = 1;
Dx1 = sqrt(1./k1).*(1i*k1.*Dx0_prop).*exp(-F*1i*k1.*rS);
Dx2 = sqrt(1./k2).*(1i*k2.*Dx0_prop).*exp(-F*1i*k2.*rS);
Dx3 = sqrt(1./k3).*(1i*k3.*Dx0_prop).*exp(-F*1i*k3.*rS);
P_field1 = zeros(size(X));
P_field2 = zeros(size(X));
P_field3 = zeros(size(X));
for xi = 1 : q : length(x0)
    progbar(1,length(x0),xi);
    R0 = sqrt((X-x0(xi)).^2 + Y.^2);
    P_field1 = P_field1 + 1/4/pi*exp(-1i*k1*R0)./R0.*Dx1(xi)*(dxSSD*q);
    P_field2 = P_field2 + 1/4/pi*exp(-1i*k2*R0)./R0.*Dx2(xi)*(dxSSD*q);
    P_field3 = P_field3 + 1/4/pi*exp(-1i*k3*R0)./R0.*Dx3(xi)*(dxSSD*q);
end
P_field{1} = P_field1;
P_field{2} = P_field2;
P_field{3} = P_field3;

y2 = (-1:dxMesh:0)';
[X2,Y2] = meshgrid(x,y2);
R0 = sqrt( (X2-xSource(1)).^2 + (Y2-xSource(2)).^2);
Pref = exp(-1i*k2*R0)./R0./4/pi;
%%
Aref1 = 1/4/pi./sqrt((X-xSource(1)).^2+(Y-xSource(2)).^2);
%
trans = 0.8;
scale = 0.2;
LW = 0.75;
LW2 = 1.5;
q= 7;
pos = [ 0.15 0.175 0.8 0.8 ];
ftsize = 13;

f = figure('Units','points','Position',[150,150,250,200]);
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(P_field2),'FaceAlpha',1);
hold on
pc2 = pcolor(x,y2,real(Pref),'FaceAlpha',0.25);
shading interp
axis equal tight
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')

cax = caxis;
hold on
fill([x0; flipud(x0)], [zeros(size(win)); win*scale], [245,154,159]/255, ...
    'FaceAlpha', trans, 'EdgeColor', 'none');
plot(x0, win*scale, 'Color',[162,50,52]/255, 'LineWidth', LW2)
hold off
clim(cax) % restore caxis


caxis([-1,1]*1/4/pi/(xRef(2)-xSource(2)))

q= 100;
ixs = find(x0<=x(end)&x0>=x(1));
ixs = ixs(1:q:end);
n0 = repmat([0,1],length(x0(ixs)));
draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.03);
off = 0.02;
line(gca().XLim,[0,0],'LineWidth',0.1,'Color',[0,0,0]);
line(gca().XLim,[1,1]*xRef(2),'LineWidth',1,'Color',[0,0,0], 'LineStyle','--');

flatIndices = find(abs(win - 1) < 1e-10);
ixEnd1 = ixs(1);
ixEnd2 = ixs(end);
ixInner1 = flatIndices(1);
ixInner2 = flatIndices(end);
xEnd1 = [x0(ixEnd1),0];
xEnd2 = [x0(ixEnd2),0];
xInner1 = [x0(ixInner1),0];
xInner2 = [x0(ixInner2),0];
k1 = [xEnd1(1)-xSource(1) 0-xSource(2)]./norm(xEnd1-xSource);
k2 = [xEnd2(1)-xSource(1) 0-xSource(2)]./norm(xEnd2-xSource);
k3 = [xInner1(1)-xSource(1) 0-xSource(2)]./norm(xInner1-xSource);
k4 = [xInner2(1)-xSource(1) 0-xSource(2)]./norm(xInner1-xSource);
D = 10;

xL= gca().XLim;
yL= gca().YLim;
line( [xSource(1), xEnd1(1)+k1(1)*D ], [xSource(2), 0+k1(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')
line( [xSource(1), xEnd2(1)+k2(1)*D ], [xSource(2), 0+k2(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')
line( [xSource(1), xInner1(1)+k3(1)*D ], [xSource(2), 0+k3(2)*D]  ,'LineWidth',LW,'Color','black','LineStyle','--')
line( [xSource(1), xInner2(1)+k4(1)*D ], [xSource(2), 0+k4(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')

xlim(xL);
ylim(yL);
xRec = [1,1.5];
xRecStat = -(xRec(1)-xSource(1))./(xRec(2)-xSource(2)).*xSource(2);
xRefStat = (xRec(1)-xSource(1))./(xRec(2)-xSource(2))*(xRef(2)-xSource(2)) + xSource(1);
scatter(xSource(1),xSource(2),20,"black",'filled')
scatter(xRec(1),xRec(2),20,"black",'filled')
scatter(xRecStat(1),0,20,"black",'filled')
scatter(xRefStat(1),xRef(2),20,"black",'filled')
line([xSource(1),xRec(1)],[xSource(2),xRec(2)],'color','black','LineWidth',1.5)

text(xRec(1)+0.075,xRec(2)+0.075, '$\mathbf{x}$','Interpreter','latex','FontSize',ftsize);
text(xRecStat(1)+0.025,0-0.175, '$\mathbf{x}_0^*(\mathbf{x})$','Interpreter','latex','FontSize',ftsize);
text(xRefStat(1)+0.075,xRef(2)-0.12, '$\mathbf{x}_{\mathrm{ref}}( \mathbf{x}_0^*(\mathbf{x}) )$','Interpreter','latex','FontSize',ftsize);
text(xSource(1)-0.25,xSource(2)+0.12, '$\mathbf{x}_{\mathrm{s}}$','Interpreter','latex','FontSize',ftsize);


set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', 'traditional_WFS_setup' ,'-dpng')
%%
% pcolor(x,y,abs(20*log10(abs(P_field)./A0)))
% shading interp
% axis equal tight
% caxis([0,3])
% colorbar

pos = [ 0.165 0.24 .8 .8];
ftsize = 13;
f = figure('Units','points','Position',[150,150,250,120]);
p2 = axes('Units','normalized','Position',pos(1,:));
X0 = -(X-xSource(1))./(Y-xSource(2)).*xSource(2);
win0 = interp1(x0,win,X0,'linear',0);
Aest = win0.*sqrt( xRef(2)./Y.*(Y-xSource(2))./(xRef(2)-xSource(2)) );
pcolor(x,y,abs(20*log10(Aest)))
shading interp
axis equal tight
clim([0,3])
cb = colorbar;
cb.Label.String = '[dB]';
hold on
cax = caxis;
hold on
fill([x0; flipud(x0)], [zeros(size(win)); win*scale], [245,154,159]/255, ...
    'FaceAlpha', trans, 'EdgeColor', 'none');
plot(x0, win*scale, 'Color',[162,50,52]/255, 'LineWidth', LW2)
hold off
clim(cax) % restore caxis

draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.03);

xL= gca().XLim;
yL= gca().YLim;
line( [xSource(1), xEnd1(1)+k1(1)*D ], [xSource(2), 0+k1(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')
line( [xSource(1), xEnd2(1)+k2(1)*D ], [xSource(2), 0+k2(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')
line( [xSource(1), xInner1(1)+k3(1)*D ], [xSource(2), 0+k3(2)*D]  ,'LineWidth',LW,'Color','black','LineStyle','--')
line( [xSource(1), xInner2(1)+k4(1)*D ], [xSource(2), 0+k4(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')

xlim(xL);
ylim(yL);
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', 'traditional_WFS_dA_est' ,'-dpng')

%%
for n = 1 : length(P_field)
    pos = [ 0.165 0.225 .8 .8];
    ftsize = 13;
    f = figure('Units','points','Position',[150,150,200,125]);
    p3 = axes('Units','normalized','Position',pos(1,:));

    pcolor(x,y,abs(20*log10(abs(P_field{n})./Aref1)))
    shading interp
    axis equal tight
    clim([0,3])
    hold on
    cax = caxis;
    hold on
    fill([x0; flipud(x0)], [zeros(size(win)); win*scale], [245,154,159]/255, ...
        'FaceAlpha', trans, 'EdgeColor', 'none');
    plot(x0, win*scale, 'Color',[162,50,52]/255, 'LineWidth', LW2)
    hold off
    clim(cax) % restore caxis

    draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.03);

    xL= gca().XLim;
    yL= gca().YLim;
    line( [xSource(1), xEnd1(1)+k1(1)*D ], [xSource(2), 0+k1(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')
    line( [xSource(1), xEnd2(1)+k2(1)*D ], [xSource(2), 0+k2(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')
    line( [xSource(1), xInner1(1)+k3(1)*D ], [xSource(2), 0+k3(2)*D]  ,'LineWidth',LW,'Color','black','LineStyle','--')
    line( [xSource(1), xInner2(1)+k4(1)*D ], [xSource(2), 0+k4(2)*D] ,'LineWidth',LW,'Color','black','LineStyle','--')

    xlim(xL);
    ylim(yL);
    xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
    ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
    set(gca,'FontName','Times New Roman');
    allAxesInFigure = findall(f,'type','axes');
    set(allAxesInFigure,'FontSize',ftsize);

    set(gcf,'PaperPositionMode','auto');
    print( '-r300', sprintf('traditional_WFS_dA_mes_%i',n) ,'-dpng')
end

%%
freq_vec = logspace(log10(10),log10(10e3),256);
taper0 = 0.5;
Lref = 1.5;
xRec0 = [(0:dxMesh:Lref)', ones(length((0:dxMesh:Lref)'),1)];

k_vec = 2*pi*freq_vec / c;
dref = 1./(abs(xRef(2)./(xRef(2)-xSource(2))))./rS;
SPA_corr = sqrt(2*pi./abs(dref)).*exp(-1i*pi/4*sign(dref));
win = tukeywin(length(x0),taper0);
Dx0_prop = 2/(4*pi)*win.*SPA_corr.*(y0-xSource(2))./rS.^2;
Dx = sqrt(1./k_vec).*(1i*k_vec.*Dx0_prop).*exp(-F*1i*k_vec.*rS);
for xi = 1 : size(xRec0,1)
    xi
    R0 = vecnorm(xRec0(xi,:)-[x0 0*x0],2,2);
    G0 = 1/4/pi*exp(-1i*k_vec.*R0 )./R0;
    P_target = 1/4/pi/norm(xRec0(xi,:)-xSource);
    P_synth(:,xi) = sum(Dx.*G0*dxSSD,1).' ./(P_target);
end
%
pos = [ 0.11 0.3 0.85 0.67];
ftsize = 13;

f = figure('Units','points','Position',[150,150,500,140]);
p1 = axes('Units','normalized','Position',pos(1,:));
semilogx(freq_vec, 20*log10(abs(P_synth(:,1))), 'k', 'LineWidth', 1.5);

hold on;

for col = 1:8:size(P_synth, 2)
    % Other columns: light gray
    semilogx(freq_vec, 20*log10(abs(P_synth(:,col))), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.05);
end

% Plot last column
sl1 = semilogx(freq_vec, 20*log10(abs(P_synth(:,1))), 'k', 'LineWidth', 2);
sl2 = semilogx(freq_vec, 20*log10(abs(P_synth(:,end))), 'k', 'LineWidth', 1);
ylim([-10,2])
grid on
xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$\frac{|P_{\mathrm{synth}}(\mathbf{x},f)|}{|P(\mathbf{x},f)|}$ [dB]', 'FontSize',ftsize,'Interpreter','latex')
legend([sl1,sl2],{'$\mathbf{x} = [0,1]$','$\mathbf{x} = [1.5,1]$'},'FontSize',ftsize,'Interpreter','latex','Location','southeast')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');

print( '-r300', sprintf('linear_wfs_transfer_taper_%i',taper0*100) ,'-dpng')
