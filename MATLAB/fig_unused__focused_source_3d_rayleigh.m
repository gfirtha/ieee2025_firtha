
clear
%close all


Lx = 10;
Lz = Lx;
dx = 0.1;
x0 = (-Lx:dx:Lx)';
z0 = (-Lz:dx:Lz)';

winX = tukeywin(length(x0),0);
winZ = tukeywin(length(z0),0);


[X0,Z0] = meshgrid(x0,z0);
Y0 = 0;
win = winZ*winX';
%win = win(:);


Nf = 250;
freq = logspace(0,log10(1e3),Nf)';
c = 343;
xRef = [0,2,0];
xS = [0,1,0];
tDirect = norm(xRef-xS)/c;
xEnd = [Lx,0,0];
tDiffr = (-norm( xEnd-xS ) + norm( xRef-xEnd ))/c;

L0 = 20*log10(norm(xRef-xS)/norm(xRef+xS));

F = sign(Y0-xS(2));
rS = sqrt( (X0-xS(1)).^2 + (Y0-xS(2)).^2 + (Z0-xS(3)).^2 );
rG = sqrt( (xRef(1)-X0).^2 + (xRef(2)-Y0).^2 + (xRef(3)-Z0).^2 );

Phi0 = exp(-1i*(F*rS+rG)/c);
A0_1 =  F*(Y0-xS(2))./rS.^2./rG;
A0_2 = (Y0-xS(2))./rS.^3./rG;

A0 = 1/(4*pi)*1/norm(xRef-xS);
%%
P1 = nan(size(freq));
P2 = nan(size(freq));
for fi = 1 : length(freq)
    fi
    freqSSDx = freq(fi);
    omega = 2*pi*freqSSDx;
    P1(fi) = -2*(1/(4*pi))^2*sum(sum( win.*(1i*omega/c*A0_1).*Phi0.^omega,1),2 ).*dx^2;
    P2(fi) = -2*(1/(4*pi))^2*sum(sum( win.*A0_2.*Phi0.^omega,1),2 ).*dx^2;
end
Pfull = P1+P2;
%%
pos = [ 0.165 0.28 .8 .7];
ftsize = 13;
f = figure('Units','points','Position',[150,150,250,150]);
p2 = axes('Units','normalized','Position',pos(1,:));

semilogx(freq,20*log10(abs(P1/A0)),'Color',[1,1,1]*0.5,'LineStyle',':','LineWidth',1.5)
hold on
semilogx(freq,20*log10(abs(P2/A0)),'Color',[1,1,1]*0.5,'LineStyle',':','LineWidth',1.5)
semilogx(freq,20*log10(abs(Pfull/A0)),'black','LineWidth',2)
xlim([freq(1),freq(end)])
grid on
ylim([-20,5])
line(gca().XLim,[1,1]*L0,'LineWidth',1.0,'Color','black','LineStyle','--');
for n = 1 : 10
    line([1,1]*(n/(tDirect-tDiffr)),gca().YLim,'LineWidth',1.0,'Color','black','LineStyle','--');
end
    

xlim([freq(1),freq(end)]);

xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$P(\mathbf{x},f)$ [dB]', 'FontSize',ftsize,'Interpreter','latex')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');
print( '-r300', sprintf( 'focused_transfer_Lx_%i_dy_%i', Lx, (xRef(2)-xS(2))) ,'-dpng')
