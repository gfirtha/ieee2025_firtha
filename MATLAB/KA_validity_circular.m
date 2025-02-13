clear
close all

c = 343;
fu = 2.5e3;
Nf = 5e3;
freq = logspace(-2,log10(fu),Nf)';
omega = 2*pi*freq;
k = omega / c;

Rssd = 2;
Nssd = 512;
fi = (0:1/Nssd:1-1/Nssd)'*2*pi;
x0 = Rssd*[cos(fi) sin(fi)];
n0 = [ -cos(fi) -sin(fi) ];
xRef = [0,0];

xSource = [3.5,0];



rhoG = sqrt(sum((xRef-x0).^2,2));
dl = Rssd*2*pi/size(x0,1);
rVS = norm(xSource);
rhoP = sqrt(sum( (x0-xSource).^2,2 ));
kh = (x0-xSource)./rhoP;
khn = sum(n0 .* kh,2);
F = -1+2*all(khn<=0); % Focused flag
if F == -1
    win = double(khn>=0);
elseif F == 1
    ks = (xSource)/norm(xSource);
    win = (kh*ks').^2.*((kh*ks')>0);
end
dref = sqrt(rhoG.*rhoP./(rhoG-F*rhoP));
tau = -F*rhoP/c;
amp = 1/(4*pi)*1./rhoP*norm(xSource)*4*pi;
tau0 =  norm(xSource)/c;
Dx = -k.*sqrt(-F*8*pi./(1i*k)).*(amp.*win.*dref.*1i.*khn.*dl)'.*exp(-1i*omega.*(tau' - tau0));

% Give analytical lower freq. limit
rP0 = min(rhoP);
rG0 = min(rhoG);
cG0 = 1./rG0;
cP0 = 1./rP0;
K1 = (2*cP0.^2+cG0.^2)./(-F.*cP0+cG0);
K2 = 3/4*(-F*cP0.^3+cG0.^3)./(-F.*cP0+cG0).^2;
omMin = c.*max( abs(K1-K2), 1./Rssd);
freqAnal =  omMin/2/pi;

% Estimate synthesized field at xRef and mesure lower freq. limit
R0 = sqrt(sum((xRef-x0).^2,2));
Ptransfer = sum( 1/(4*pi)*Dx.*exp(-1i*omega*R0'/c)./R0',2);
L0 = -1;
ixs = find(20*log10(abs(Ptransfer))<L0+0.1 & 20*log10(abs(Ptransfer))>L0-0.1);
if ~isempty(ixs)
    freqMeasured = freq(ixs(1));
else
    freqMeasured = 0;
end

%%
dx = 0.05;
f0 = 1000;
[~,ix] = min(abs(freq-f0));
f0 = freq(ix);
Dx0 = Dx(ix,:);
x = (-3:dx:3)';
y = (-3:dx:3)';
[X,Y] = meshgrid(x,y);
Psynth = zeros(size(X));
%
k0 = 2*pi*f0/c;
for n = 1 : length(Dx0)
    n
    R0 = sqrt( (X-x0(n,1)).^2 + (Y-x0(n,2)).^2 );
    Psynth = Psynth + 1/4/pi*exp(-1i*k0*R0)./R0.*Dx0(n);
end
%%
figure
subplot(2,1,1)
pcolor(x,y,real(Psynth))
shading interp
caxis([-1,1])
axis equal tight

hold on
scatter(xSource(1),xSource(2),20,'red','filled');
scatter(xRef(1),xRef(2),20,'black','filled');
grid on
axis equal tight
xlim([x(1),x(end)])

subplot(2,1,2)
semilogx(freq, 20*log10(abs(Ptransfer)))
grid on
hold on
xlim([freq(1),freq(end)])
ylim([-20,20 ])
yl = get(gca,'ylim');
hf1 = plot([freqAnal,freqAnal],yl,'--r','LineWidth',1.25);
hf2 = plot([freqMeasured,freqMeasured],yl,'--g','LineWidth',1.25);
legend([hf1, hf2], {'f_{analytical}', 'f_{measured}'})
xlabel('f [Hz]')
ylabel('P_{synth} [dB]')