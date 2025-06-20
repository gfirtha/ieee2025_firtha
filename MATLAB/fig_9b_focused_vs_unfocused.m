clear
%close all

rSSD = 2;         % SSD radius
nSSD = 360*3;       % secondary source number
xRef = [0,0];
c = 343;
taper = 0.25;
freq = logspace(log10(10),log10(20e3),256);
k = 2*pi*freq/c;
phi = (0:nSSD-1)'/nSSD'*2*pi + pi;
dphi = mean(diff(phi));
x0 = rSSD*[cos(phi) sin(phi)];
n0 = [ -cos(phi) -sin(phi) ];
rhoG = sqrt(sum(x0.^2,2));
dl = mean(rhoG*2*pi/size(x0,1));

xSource = [ 1, 0 ];    % source position
rhoP = sqrt(sum( (x0-xSource).^2,2 ));
kh = (x0-xSource)./rhoP;
khn = sum(n0 .* kh,2);
F = -(-1+2*all(khn<0)); % Focused flag
if F == 1
    win = double(khn>=0);
elseif F == -1
    ks = (xSource)/norm(xSource);
    %    win = (kh*ks').^2.*((kh*ks')>0);
    win = (kh*ks').^0.*((kh*ks')>0);
end
[~,max_ix] = max(abs(khn.*win));

for n = 1 : 4
    win0 = tukeywin(length(find(win==1)),0.25*n);
    win0 = [win0;zeros(length(find(win~=1)),1)];
    win0 = circshift(win0,max_ix-round(length(find(win==1))/2)-1);
    win1{n} = win.*win0;

    rRef = rhoG;
    dref = sqrt(rRef.*rhoP./(rRef+F*rhoP));
    Px0 = 1/(4*pi)*exp(-1i*F*k.*rhoP)./rhoP;
    Dx = sqrt(F*8*pi*1i*k).*(win1{n}.*dref.*khn.*Px0*dl);
    Psynth{n} = sum(Dx.*exp(-1i*k*rSSD)./rSSD/4/pi,1) / (1/4/pi/norm(xSource));
end

alpha = linspace(0.25,1,16);
for n = 1 : length(alpha)
    win0 = tukeywin(length(find(win==1)),alpha(n));
    win0 = [win0;zeros(length(find(win~=1)),1)];
    win0 = circshift(win0,max_ix-round(length(find(win==1))/2)-1);
    win1{n} = win.*win0;

    rRef = rhoG;
    dref = sqrt(rRef.*rhoP./(rRef+F*rhoP));
    Px0 = 1/(4*pi)*exp(-1i*F*k.*rhoP)./rhoP;
    Dx = sqrt(F*8*pi*1i*k).*(win1{n}.*dref.*khn.*Px0*dl);
    Psynth2{n} = sum(Dx.*exp(-1i*k*rSSD)./rSSD/4/pi,1) / (1/4/pi/norm(xSource));
end

xSource = [ 3, 0 ];    % source position
rhoP = sqrt(sum( (x0-xSource).^2,2 ));
kh = (x0-xSource)./rhoP;
khn = sum(n0 .* kh,2);
win = double(khn>=0);
dref = sqrt(rhoG.*rhoP./(rhoG+rhoP));
Px0 = 1/(4*pi)*exp(-1i*k.*rhoP)./rhoP;
Dx = sqrt(8*pi*1i*k).*(win.*dref.*khn.*Px0*dl);
Psynth{5} = sum(Dx.*exp(-1i*k*rSSD)./rSSD/4/pi,1) / (1/4/pi/norm(xSource));

ftsize = 13;
LW = 1.5;
colors = linspace(0.25,0.9,6)'*[1 1 1];

pos = [ 0.11 0.3 0.85 0.67];
ftsize = 13;

f = figure('Units','points','Position',[150,150,500,140]);
p1 = axes('Units','normalized','Position',pos(1,:));
semilogx(freq,20*log10(abs(Psynth{1})), 'Color',colors(2,:),'LineWidth',LW)
hold on
for n = 1 : length(Psynth2)
    semilogx(freq,20*log10(abs(Psynth2{n})), 'Color',[1,1,1]*0.8,'LineWidth',0.5)

end
xlim([freq(1), freq(end)])
s1 = semilogx(freq,20*log10(abs(Psynth{1})), 'Color',colors(2,:),'LineWidth',LW)
s2 = semilogx(freq,20*log10(abs(Psynth{2})), 'Color',colors(3,:),'LineWidth',LW)
s3 = semilogx(freq,20*log10(abs(Psynth{3})), 'Color',colors(4,:),'LineWidth',LW)
s4 = semilogx(freq,20*log10(abs(Psynth{4})), 'Color',colors(5,:),'LineWidth',LW)
s5 = semilogx(freq,20*log10(abs(Psynth{5})),'Color',colors(1,:),'LineWidth',LW,'LineStyle','--')
grid on

xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$\frac{|P_{\mathrm{synth}}(\mathbf{0},f)|}{|P(\mathbf{0},f)|}$ [dB]', 'FontSize',ftsize,'Interpreter','latex')
legend([s1,s2,s3,s4,s5],{'foc., $\alpha = 0.25$',...
    'foc., $\alpha = 0.5$',...
    'foc., $\alpha = 0.75$',...
    'foc., $\alpha = 1$','unfocused'}...
    ,'FontSize',ftsize-3,'Interpreter','latex','Location','southeast')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'focused_vs_unfocused_freq','-dpng')

function r = dist(x)

r = sqrt( sum( x.^2 , 2)  );
end
