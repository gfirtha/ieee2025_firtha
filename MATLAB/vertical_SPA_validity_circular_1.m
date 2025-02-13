clear
close all
R0 = 2;
fi = linspace(-pi/4,pi/4,91)'-pi/2;

x0_vec = R0*[cos(fi) sin(fi)];
dx = 0.05;
z0 = (-100:dx:100)';


xSource = [0,-3,0];
xRef = [0,0,0];

fu = 5e3;
Nf = 2e3;
freq = logspace(-1,log10(fu),Nf)';

c = 343;
omega = 2*pi*freq;
k = omega/c;

F = sign(-R0 - xSource(2));
for xi = 1 : length(x0_vec)
    xi
    x0 = x0_vec(xi,1);
    y0 = x0_vec(xi,2);


    % Numerical integral
    rP = sqrt((x0-xSource(1)).^2+(y0-xSource(2)).^2 + (z0-xSource(3)).^2);
    rG = sqrt((x0-xRef(1)).^2+(y0-xRef(2)).^2 + (z0-xRef(3)).^2);
    A = 1./rP.^2./rG;
    Phi = -(F*rP+rG);
    I = sum( A.'.*exp(1i*k*Phi.'), 2) * dx;

    % SPA
    rP0 = sqrt((x0-xSource(1)).^2+(y0-xSource(2)).^2);
    rG0 = sqrt((x0-xRef(1)).^2+(y0-xRef(2)).^2);
    A0 = 1./rP0.^2./rG0;
    Phi0 = -(F*rP0+rG0);
    ddPhi0 = (F./rP0 + 1./rG0);
    Iappr = sqrt(2*pi./(1i*k*abs(ddPhi0)))*exp(1i*pi/4*sign(ddPhi0))*A0.'.*exp(1i*k*Phi0.');


    cG0 = 1./rG0;
    cP0 = 1./rP0;
    K1 = (2*cP0.^2+cG0.^2)./(F.*cP0+cG0);
    K2 = 3/4*(F*cP0.^3+cG0.^3)./(F*cP0+cG0).^2;

    omMin = c.*abs(K1-K2);
    freqAnal(xi) = omMin / 2 / pi;

    dBdiff = abs(20*log10(abs(I))-20*log10(abs(Iappr)));
    ixs = find(dBdiff<1.5);
    freqMeasured(xi) = freq(ixs(1)) ;

end
%%
figure
plot(fi,freqMeasured)
hold on
plot(fi,freqAnal,'--')
legend('Measured cutott','Analytical cutoff')
xlabel('\phi [rad]')
ylabel('f [Hz]')
xlim([fi(1),fi(end)])
grid on