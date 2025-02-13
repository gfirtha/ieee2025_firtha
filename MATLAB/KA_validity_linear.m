clear
close all
Lx = 50;
dx = 0.05;

x0 = (-Lx:dx:Lx);
y0 = 0;
fu = 2.5e3;
Nf = 5e3;
freq = logspace(-2,log10(fu),Nf)';

c = 343;
omega = 2*pi*freq;
k = omega/c;

xSource = [0,-1];
xRef = [0,2];

F = -sign(xSource(2));

rP = sqrt( (x0-xSource(1)).^2 + (y0-(xSource(2))).^2 );
rG = sqrt( (x0-xRef(1)).^2 + (y0-xRef(2)).^2 );

% Linear reference
dref = abs(xRef(2)./(xRef(2)-xSource(2)));
Dx = -1/4/pi*sqrt(8*pi/1i./k).*sqrt(dref).*1i.*k.*xSource(2)./sqrt(rP).*exp(-F*1i*k.*rP)./rP;

% 2.5D Kirchhof
% kn = -F*k*(y0-xSource(2))./rP;
% P0 = 1/4/pi*exp(-F*1i*k.*rP)./rP;
% dref = 1./(abs(1./rP + 1./rG));
% Dx = 2*sqrt(2*pi./(1i*k).*dref).*1i.*kn.*P0;
% 
%
% Bleistein p. 130

P = zeros(length(freq),1);
R0 = sqrt( (x0-xRef(1)).^2 + (y0-xRef(2)).^2 );
for n = 1 : size(Dx,2)
    P = P + 1/4/pi*dx*exp(-1i*k*R0(n))./R0(n).*Dx(:,n)*norm(xRef-xSource)*4*pi;
end
%

xi = line_intersection(xSource', xRef', [x0(1),y0]', [x0(end),y0]')';
[~, stat_ix] = min(abs(xi(1)-x0));

rP0 = sqrt( (xi(1)-xSource(1)).^2 + (xi(2)-(xSource(2))).^2 );
rG0 = sqrt( (xi(1)-xRef(1)).^2 + (xi(2)-xRef(2)).^2 );
cG0 = 1./rG0;
cP0 = 1./rP0;
K1 = (2*cP0.^2+cG0.^2)./(F.*cP0+cG0);
K2 = 3/4*(F*cP0.^3+cG0.^3)./(F.*cP0+cG0).^2;

omMin = c.*abs(K1-K2);
freqAnal = omMin / 2 / pi;

%
f0 = 1000;
[~,ix] = min(abs(freq-f0));
f0 = freq(ix);
Dx0 = Dx(ix,:);
x = (-3:dx:3)';
y = (0:dx:4)';
[X,Y] = meshgrid(x,y);

Psynth = zeros(size(X));
%
k = 2*pi*f0/c;
for n = 1 : length(Dx0)
    n
    R0 = sqrt( (X-x0(n)).^2 + Y.^2 );
    Psynth = Psynth + 1/4/pi*dx*exp(-1i*k*R0)./R0.*Dx0(n)*norm(xRef-xSource)*4*pi;
end
%%
xData = log10(freq);
yData = 20*log10(abs(P));

L0 = -1;
ixs = find(20*log10(abs(P))<L0+0.1 & 20*log10(abs(P))>L0-0.1);
if ~isempty(ixs)
    freqMeasured = freq(ixs(1));
else
    freqMeasured = 0;
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
scatter(xi(1),xi(2),20,'black','filled');
grid on
axis equal tight
line( [xSource(1),xRef(1)], [xSource(2),xRef(2)],'Color','black')
line( [x0(1),x0(end)], [y0,y0],'Color','black')
xlim([x(1),x(end)])

subplot(2,1,2)
semilogx(freq, 20*log10(abs(P)))
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



%%
