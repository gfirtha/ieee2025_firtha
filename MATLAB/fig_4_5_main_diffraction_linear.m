clear
%close all
addpath(genpath( 'Functions\' ))

Lssd = 1.5;                 % Half-length of integration domain in x-direction
dx = 1e-2;
taper = 0.5;
freq_vec = logspace(log10(100),log10(10e3),256); % Frequency vector (logarithmic scale)
c = 343;                    % Speed of sound in air (m/s)
k_vec = 2*pi*freq_vec / c;
xRef = [0,1];       % Receiver/reference position
xSource = [0,-1e3];        % Source position

x0 = (-20:dx:20)';
y0 = 0;                   % Fixed y-position of the integration line

ixs = find(x0>=-Lssd & x0<=Lssd);
win = zeros(length(x0),1);
win(ixs) = tukeywin(length(ixs),taper);

win_ = (1-win).*tukeywin(length(win),0.2);
P_0 = 1/4/pi/norm(xSource-xRef);

Lmes = 5;
xMes = (-Lmes:dx:Lmes)';
xRec = [xMes, xRef(2)*ones(length((-Lmes:dx:Lmes)),1)];

f0 = 1e3;
k0 = 2*pi*f0 / c;
Ly = 3;
y = (0:dx:Ly);
[X,Y] = meshgrid(xMes,y);

[~, fix] = min(abs(freq_vec-f0));
[~,xix] = min(abs(xRef(2)-y));
Lx = Lmes;


F = sign(y0-xSource(2));
rS = sqrt( (x0-xSource(1)).^2 + (y0-xSource(2)).^2 );

SPA_corr = sqrt(abs(xRef(2)./(xRef(2)-xSource(2))) );
Dx0_prop = 1/sqrt(2*pi)*SPA_corr.*(y0-xSource(2));
Dx = sqrt(1i*k_vec).*Dx0_prop.*exp(-F*1i*k_vec.*rS)./rS.^1.5;
Dx_f0 = sqrt(1i*k0).*Dx0_prop.*exp(-F*1i*k0.*rS)./rS.^1.5;

x0_stat = -(X-xSource(1))./(Y-xSource(2)).*xSource(2);
flatIndices = find(abs(win - 1) < 1e-10);
ixEnd1 = ixs(1);
ixEnd2 = ixs(end);
ixInner1 = flatIndices(1);
ixInner2 = flatIndices(end);
xEnd1 = [x0(ixEnd1),0];
xEnd2 = [x0(ixEnd2),0];
xInner1 = [x0(ixInner1),0];
xInner2 = [x0(ixInner2),0];

domain_shadow = (x0_stat<-Lssd) | (x0_stat>Lssd);
domain_inter = ((x0_stat<=x0(ixInner1) & x0_stat>-Lssd) ) | ((x0_stat>=x0(ixInner2) & x0_stat<Lssd));
domain_illum    = ( x0_stat>=x0(ixInner1) & x0_stat<=x0(ixInner2));

Iw0 =  interp1(x0,win,x0_stat,'linear',0);
Aest = Iw0.*sqrt( xRef(2)./Y.*(Y-xSource(2))./(xRef(2)-xSource(2)) );

R = sqrt((X-xSource(1)).^2 + (Y-xSource(2)).^2);
P_target_xw = Aest(xix,:).'./4/pi./vecnorm(xRec-xSource,2,2).*exp(-1i*vecnorm(xRec-xSource,2,2).*k_vec);
P_target_f0 = Aest./4/pi.*exp(-1i*k0*R)./R;

P_synth_xw = zeros(length(xRec),length(freq_vec));
IA_xw = zeros(length(xRec),length(freq_vec));
P_synth_f0 = zeros(size(X));
IA = zeros(size(X));
for xi = 1 : length(x0)
    xi
    R0 = vecnorm(xRec-[x0(xi,1) 0],2,2);
    G0 = 1/4/pi*exp(-1i*k_vec.*R0 )./R0;
    P_synth_xw = P_synth_xw + win(xi).*Dx(xi,:).*G0*dx;
    IA_xw = IA_xw + win_(xi).*Dx(xi,:).*G0*dx;

    R0 = sqrt( (X-x0(xi)).^2 + (Y-y0).^2 );
    G0 = 1/4/pi*exp(-1i*k0.*R0 )./R0;
    P_synth_f0 = P_synth_f0 + win(xi).*Dx_f0(xi).*G0*dx;
    IA = IA + win_(xi).*Dx_f0(xi).*G0*dx;

end
P_diffr_num_xw =  (P_synth_xw-P_target_xw).*domain_inter(xix,:).'-...
    IA_xw.*domain_illum(xix,:).' + ...
    P_synth_xw.*domain_shadow(xix,:).';

%P_diffr_num_xw =  (P_synth_xw-P_target_xw).*(domain_inter(xix,:).' - domain_illum(xix,:).') + P_synth_xw.*domain_shadow(xix,:).';
P_diffr_f0 =  (P_synth_f0-P_target_f0).*domain_inter - IA.*domain_illum + domain_shadow.*P_synth_f0;
%%
trans = 0.8;
scale = 0.5;
LW = 0.5;
LW2 = 1.5;
q= 16;

P_comp = {P_synth_f0,P_target_f0,P_diffr_f0};
pos = [ 0.125 0.2 0.85 0.8 ];
ftsize = 13;
for n = 1 : length(P_comp)
    f = figure('Units','points','Position',[150,150,250,150]);
    p1 = axes('Units','normalized','Position',pos(1,:));
    pcolor(X,Y,(real(P_comp{n}/P_0)));
    axis equal tight
    shading flat
    xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
    ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
    xlim([-1,1]*3)
    caxis([-1,1]*0.5)
    cax = caxis;
    hold on
    fill([x0; flipud(x0)], [zeros(size(win)); win*max(cax)*scale], [245,154,159]/255, ...
        'FaceAlpha', trans, 'EdgeColor', 'none');
    plot(x0, win*max(cax)*scale, 'Color',[162,50,52]/255, 'LineWidth', LW2)
    hold off
    clim(cax) % restore caxis
    yl = get(gca, 'YLim');
    xl = get(gca, 'XLim');
    line([xEnd1(1),xEnd1(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
    line([xEnd2(1),xEnd2(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
    line([xInner1(1),xInner1(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
    line([xInner2(1),xInner2(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
    line([xl(1),xEnd1(1)],[0,0],'Color','black','LineStyle','-','LineWidth',LW2)
    line([xEnd2(1),xl(2)],[0,0],'Color','black','LineStyle','-','LineWidth',LW2)
    ixs = find(x0<=xEnd2(1)&x0>=xEnd1(1));
    ixs = ixs(1:q:end);
    n0 = repmat([0,1],length(x0(ixs)));
    draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.05);
    set(gca,'FontName','Times New Roman');
    allAxesInFigure = findall(f,'type','axes');
    set(allAxesInFigure,'FontSize',ftsize);
    set(gcf,'PaperPositionMode','auto');
    print( '-r300', sprintf('P_components_re_%i',n) ,'-dpng');
end
%%
pos = [ 0.1 0.2 0.835 0.8 ];
ftsize = 13;
for n = 1 : length(P_comp)

f = figure('Units','points','Position',[150,150,310,150]);

p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(X,Y,20*log10(abs(P_comp{n}/P_0)));
axis equal tight
shading flat
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')

clim([-50,1])
xlim([-1,1]*3)
cax = caxis;
hold on
scale2 = scale/2;
fill([x0; flipud(x0)], [zeros(size(win)); win*scale2], [245,154,159]/255, ...
    'FaceAlpha', trans, 'EdgeColor', 'none');
plot(x0, win*scale2, 'Color',[162,50,52]/255, 'LineWidth', LW2);
hold off
caxis(cax) % restore caxis
line([xEnd1(1),xEnd1(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xEnd2(1),xEnd2(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xInner1(1),xInner1(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xInner2(1),xInner2(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xl(1),xEnd1(1)],[0,0],'Color','black','LineStyle','-','LineWidth',LW2)
line([xEnd2(1),xl(2)],[0,0],'Color','black','LineStyle','-','LineWidth',LW2)
cb = colorbar;
cb.Label.String = '[dB]';
ixs = find(x0<=xEnd2(1)&x0>=xEnd1(1));
ixs = ixs(1:q:end);
n0 = repmat([0,1],length(x0(ixs)));
draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.05);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
    print( '-r300', sprintf('P_components_abs_%i',n) ,'-dpng');
end

%%
ixCritical = [ixEnd1(1), ixEnd2(1), ixInner1(1), ixInner2(1)]';
xCritical = x0(ixCritical);
P_diffr_anal_xw = zeros(length(xRec),length(freq_vec));
P_diffr_anal_mesh = zeros(size(X));

for i = 1 : length(xCritical) 
    Rre_ref = vecnorm([xCritical(i),0]-xRec,2,2);
    dkx = (xCritical(i)-xRec(:,1))./Rre_ref + (xCritical(i)-xSource(1))./norm( [xCritical(i),0]-xSource);
    a = (-1)^i*2*pi^2/taper^2/(2*Lssd)^2;
    P_diffr_anal_xw = P_diffr_anal_xw +  (1./(1i*k_vec.*dkx)).^3.*a.*Dx(ixCritical(i),:).*exp(-1i*k_vec.*Rre_ref)./Rre_ref./4/pi;

    Rre_mesh = sqrt( (X-xCritical(i)).^2 + Y.^2 );
    dkx_mesh = (xCritical(i)-X)./Rre_mesh + (xCritical(i)-xSource(1))./norm( [xCritical(i),0]-xSource);
    P_diffr_anal_mesh = P_diffr_anal_mesh +  (1./(1i*k0.*dkx_mesh)).^3.*a.*Dx(ixCritical(i),fix).*exp(-1i*k0.*Rre_mesh)./Rre_mesh./4/pi;

end    
%%
[~,fix2] = min(abs(freq_vec - 5e3));
[~,ixCenter] = min(abs(xRec(:,1)-0));

%%
pos = [ 0.11 0.2 0.835 0.8 ];
ftsize = 13;
f = figure('Units','points','Position',[150,150,310,150]);
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(X,Y,20*log10(abs(P_diffr_anal_mesh/P_0)));
axis equal tight
shading flat
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
yl = get(gca, 'YLim');
xl = get(gca, 'XLim');
clim([-50,1])
xlim([-1,1]*3)
cax = caxis;
hold on
scale2 = scale/2;
fill([x0; flipud(x0)], [zeros(size(win)); win*scale2], [245,154,159]/255, ...
    'FaceAlpha', trans, 'EdgeColor', 'none');
plot(x0, win*scale2, 'Color',[162,50,52]/255, 'LineWidth', LW2);
hold off
caxis(cax) % restore caxis
line([xEnd1(1),xEnd1(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xEnd2(1),xEnd2(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xInner1(1),xInner1(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xInner2(1),xInner2(1)],[yl(1),yl(2)],'Color','black','LineStyle','--','LineWidth',LW)
line([xl(1),xEnd1(1)],[0,0],'Color','black','LineStyle','-','LineWidth',LW2)
line([xEnd2(1),xl(2)],[0,0],'Color','black','LineStyle','-','LineWidth',LW2)
cb = colorbar;
cb.Label.String = '[dB]';
ixs = find(x0<=xEnd2(1)&x0>=xEnd1(1));
q = 16;
ixs = ixs(1:q:end);
n0 = repmat([0,1],length(x0(ixs)));
draw_ssd( [], [x0(ixs),0*x0(ixs)] , n0, 0.05);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'P_diffr_tukey' ,'-dpng');

%%
pos = [ 0.155 0.25 0.725 0.75];
ftsize = 13;
f = figure('Units','points','Position',[150,150,310,160]);

p1 = axes('Units','normalized','Position',pos(1,:));
semilogx(freq_vec, 20*log10(abs(P_diffr_num_xw(ixCenter,:))/P_0),'-k', 'LineWidth', 1.5);
grid on
hold on
semilogx(freq_vec, 20*log10(abs(P_diffr_anal_xw(ixCenter,:))/P_0),':','Color',[1,1,1]*0.5, 'LineWidth', 2);
xlim([freq_vec(1) freq_vec(end)])

xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$|P_{\mathrm{diffr}}(0,y_{\mathrm{ref}})|$ [dB]', 'FontSize',ftsize,'Interpreter','latex')
legend('Numerical result','Analytical model','Interpreter','latex')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'P_diffr_f_tukey' ,'-dpng');
%%
freqs_to_plot = [800, 2000, 5000]; % in Hz
shifts_dB = [50, 0, -50];          % Vertical offsets
ix_plot = arrayfun(@(f) find(abs(freq_vec - f) == min(abs(freq_vec - f)), 1), freqs_to_plot);

colors = lines(length(freqs_to_plot));
colors = [0 0.45 0.6]'*[1 1 1]; % override first color (black)

pos = [ 0.08 0.12 0.9 0.87];
ftsize = 13;
f = figure('Units','points','Position',[150,150,270,310]);

p1 = axes('Units','normalized','Position',pos(1,:));
hold on

xL = [-3, 3];
yL = [-170, 100];
for i = 1:length(ix_plot)
    ix = ix_plot(i);
    shift = shifts_dB(i);

    dB_num = 20*log10(abs(P_diffr_num_xw(:,ix) / P_0)) + shift;
    dB_ana = 20*log10(abs(P_diffr_anal_xw(:,ix) / P_0)) + shift;

    plot(xRec(:,1), dB_num, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(xRec(:,1), dB_ana, ':', 'Color', colors(i,:), 'LineWidth', 1.5);

    % Left label (true dB)
    [~,ixLeft] = min(abs(xRec(:,1)-xL(1)));
    true_dB = 20*log10(abs(P_diffr_num_xw(ixLeft,ix) / P_0));
    text(xL(1), dB_num(ixLeft)+20, sprintf('%.0f dB', true_dB), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', 10, 'Color', colors(i,:),'Interpreter','latex');

    % Right label (frequency)
    [~,ixRight] = min(abs(xRec(:,1)-xL(2)));
    text(xL(2), dB_num(ixRight)-8, sprintf('%d Hz', freqs_to_plot(i)), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'Color', colors(i,:),'Interpreter','latex');
end
legend('Numerical result','Analytical model','location','northeast','FontSize',ftsize-2,'Interpreter','latex')
% Axes settings
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$|P_{\mathrm{diffr}}(x,y_{\mathrm{ref}},f)|$ [dB]', 'Interpreter','latex');
set(gca,'FontSize',12);
xlim(xL)
ylim(yL)
grid on
set(gca, 'YTickLabel', []) % hide labels but keep ticks and grid

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'P_diffr_xf_tukey' ,'-dpng');
%%

