clear
%close all
dx = 0.01;

x0 = 0;
y0 = 0;
z0 = (-100:dx:100)';


xSource = [0,-1,0];
xRef = [0,1,0];

fu = 10e3;
Nf = 5e3;
freq = logspace(0,log10(fu),Nf)';

c = 343;
omega = 2*pi*freq;
k = omega/c;

F = -sign(xSource(2));

% Numerical integral
rP = sqrt((x0-xSource(1)).^2+(y0-xSource(2)).^2 + (z0-xSource(3)).^2);
rG = sqrt((x0-xRef(1)).^2+(y0-xRef(2)).^2 + (z0-xRef(3)).^2);
A = 1./rP.^2./rG;
Phi = -(F*rP+rG);
I = sum( A.'.*exp(1i*k*Phi.'), 2) * dx;

% SPA
rP0 = sqrt((x0-xSource(1)).^2+(y0-xSource(2)).^2 + (0-xSource(3)).^2);
rG0 = sqrt((x0-xRef(1)).^2+(y0-xRef(2)).^2 + (0-xRef(3)).^2);
A0 = 1./rP0.^2./rG0;
Phi0 = -(F*rP0+rG0);
ddPhi0 = (F./rP0 + 1./rG0);
Iappr = sqrt(2*pi./(1i*k*ddPhi0))*exp(1i*pi/4*sign(ddPhi0))*A0.'.*exp(1i*k*Phi0.');


cG0 = 1./rG0;
cP0 = 1./rP0;
K1 = (2*cP0.^2+cG0.^2)./(F.*cP0+cG0);
K2 = 3/4*(F*cP0.^3+cG0.^3)./(F.*cP0+cG0).^2;

omMin = c.*abs(K1-K2)/24;
freqAnal  = omMin / 2 / pi;

dBdiff = abs(20*log10(abs(I))-20*log10(abs(Iappr)));
ixs = find(dBdiff< 1.5);
freqMeasured = freq(ixs(1)) ;

figure
semilogx(freq, 20*log10(abs(I.*sqrt(omega))),'LineWidth',1.25);
hold on
semilogx(freq, 20*log10(abs(Iappr.*sqrt(omega))),'LineWidth',1.25);
grid on
xlim([freq(1),freq(end)])
yl = get(gca,'ylim');
hf2 = plot([freqMeasured,freqMeasured],yl,'-k','LineWidth',1.25);
hf1 = plot([freqAnal,freqAnal],yl,'--g','LineWidth',1.25);
legend('Numerical integral', 'SPA result', 'Measured frequency limit','Analytical frequency limit','location','best')
xlabel('$f$ [Hz]','Interpreter','latex')
ylabel('$20 \log10( |I| )$ [dB]','Interpreter','latex')


%%
% xl = [-1.5, 1.5];
% yl = [-1.2, 1.2];
% 
% dx = 0.01;
% x = (xl(1):dx:xl(end))';
% y = (yl(1):dx:yl(end))';
% [X,Y] = meshgrid(x,y);
% R = sqrt((X-xSource(1)).^2 + (Y-xSource(2)).^2);
% P = exp(-1i*2*pi*1e3/c*R)./R;
% 
% alphaMap = 0.5 * ones(size(P));    % start with alpha=0.5 everywhere
% alphaMap(Y > 0) = 0.5;            % more transparent above y=0
% alphaMap(Y < 0) = 0.2;            % less transparent below y=0
% 
% figure;
% p = pcolor(x,y,real(P));
% colormap gray;
% caxis([-1,1]);
% 
% % Make sure shading is flat (per face)
% shading flat
% 
% % Now you can set face-based alpha
% p.FaceColor = 'flat';
% p.FaceAlpha = 'flat';
% p.AlphaData = alphaMap;          % your 2D alpha matrix
% p.AlphaDataMapping = 'none';
% 
% hold on
% scatter(xSource(1),xSource(2),30,'red','filled');
% text(xSource(1)+0.1,xSource(2) - 0.15,'virtual source position','Interpreter','latex')
% scatter(xRef(1),xRef(2),30,'black','filled');
% text(xRef(1)+0.1,xRef(2) - 0.1,'receiver position','Interpreter','latex')
% 
% axis equal tight
% grid on
% %line( [xSource(1),xRef(1)], [xSource(2),xRef(2)],'Color','black')
% line( [x0(1)-2,x0(end)+2], [y0,y0],'Color','black','lineWidth',2 );
% text(x(1)+0.1,0.1,'SSD','Interpreter','latex')
% xlim(xl)
% ylim(yl)
% xlabel('$x$ [m]','Interpreter','latex')
% ylabel('$y$ [m]','Interpreter','latex')
% 
% %view(3)