clear
close all

c = 343;
fu = 2.5e3;
Nf = 5e3;
freq = logspace(-2,log10(fu),Nf)';
omega = 2*pi*freq;
k = omega / c;

R0 = 2;
Nssd = 512;
fi = (0:1/Nssd:1-1/Nssd)'*2*pi;
x0 = R0*[cos(fi) sin(fi)];
n0 = [ -cos(fi) -sin(fi) ];
xRef = [0,0];

rhoG = sqrt(sum((xRef-x0).^2,2));
dl = R0*2*pi/size(x0,1);

dx2 = 0.1;
xs = (0:dx2:4)';
ys = (0:dx2:3)';
[xSource,ySource] = meshgrid(xs,ys);
xSource = [xSource(:),ySource(:)];
%%
for si = 1 : size(xSource)
    si
    xSourceAct = xSource(si,:);

    rVS = norm(xSourceAct);
    rhoP = sqrt(sum( (x0-xSourceAct).^2,2 ));
    kh = (x0-xSourceAct)./rhoP;
    khn = sum(n0 .* kh,2);
    F = -1+2*all(khn<=0); % Focused flag
    if F == -1
        win = double(khn>=0);
    elseif F == 1
        ks = (xSourceAct)/norm(xSourceAct);
        win = (kh*ks').^2.*((kh*ks')>0);
    end
    dref = sqrt(rhoG.*rhoP./(rhoG-F*rhoP));
    tau = -F*rhoP/c;
    amp = 1/(4*pi)*1./rhoP*norm(xSourceAct)*4*pi;
    tau0 =  norm(xSourceAct)/c;
    Dx = -k.*sqrt(-F*8*pi./(1i*k)).*(amp.*win.*dref.*1i.*khn.*dl)'.*exp(-1i*omega.*(tau' - tau0));

    % Give analytical lower freq. limit
    rP0 = min(rhoP);
    rG0 = min(rhoG);
    cG0 = 1./rG0;
    cP0 = 1./rP0;
    K1 = (2*cP0.^2+cG0.^2)./(-F.*cP0+cG0);
    K2 = 3/4*(-F*cP0.^3+cG0.^3)./(-F.*cP0+cG0).^2;
    omMin = c.*abs(K1-K2);

    rS = norm(xSourceAct);
    K = 1/4/F*abs( abs(R0-rS)./R0./rS - 5./(R0-rS) );
    omMin = c.*max( K, 2./R0);

    freqAnal(si) =  omMin/2/pi;

    % Estimate synthesized field at xRef and mesure lower freq. limit
    Rssd = sqrt(sum((xRef-x0).^2,2));
    Ptransfer = sum( 1/(4*pi)*Dx.*exp(-1i*omega*Rssd'/c)./Rssd',2);
    L0 = -1;
    ixs = find(20*log10(abs(Ptransfer))<L0+0.1 & 20*log10(abs(Ptransfer))>L0-0.1);
    if ~isempty(ixs)
        freqMeasured(si) = freq(ixs(1));
    else
        freqMeasured(si) = 0;
    end

end
%%
%%
freqMapMeasured = reshape(freqMeasured,size(ys,1),size(xs,1));
freqMapMeasuredTemp = freqMapMeasured;
freqMapMeasuredTemp(freqMapMeasuredTemp<1e-9) = Inf;
freqMapAnal = reshape(freqAnal,size(ys,1),size(xs,1));


[xOut,yOut, freqMapAnalFull] = createFullMx(xs,ys,freqMapAnal);
[xOut,yOut, freqMapMeasuredFull] = createFullMx(xs,ys,freqMapMeasuredTemp);

figure
subplot(1,3,1)
s1 = surf(xOut,yOut,freqMapMeasuredFull);
xlabel('x_s [m]')
ylabel('y_s [m]')
title('Measured low frequency limit')
zlim([0,600])
caxis([0,1000])

subplot(1,3,2)
s2 = surf(xOut,yOut,(freqMapAnalFull));
xlabel('x_s [m]')
ylabel('y_s [m]')
title(['Analyitical low frequency limit'])
zlim([0,600])
caxis([0,1000])

subplot(1,3,3)
plot(ys,freqMapMeasured(:,1))
hold on
plot(ys,freqMapAnal(:,1))
xlabel('x_s [m]')
ylabel('f_{LF} [Hz]')
grid on
%ylim([0,600])
title(['Low frequency limit over x_s = 0'])
legend('f_{measured}',['f_{analytical} '])

