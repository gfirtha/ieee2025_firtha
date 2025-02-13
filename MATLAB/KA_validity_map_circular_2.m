clear
close all

c = 343;
fu = 2.5e3;
Nf = 10e3;
freq = logspace(-2,log10(fu),Nf)';
omega = 2*pi*freq;
k = omega / c;

Rssd = 2;
Nssd = 512;
fi = (0:1/Nssd:1-1/Nssd)'*pi-pi/2;
x0 = Rssd*[cos(fi) sin(fi)];
n0 = [ -cos(fi) -sin(fi) ];
xRef = [0,0];

rhoG = sqrt(sum((xRef-x0).^2,2));
dl = Rssd*2*pi/size(x0,1)/2;

dx2 = 0.05;
xs = (0:dx2:4)';
%xs = 1.55;
ys = 0;
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
    F = 1-2*all(khn<=0); % Focused flag
    if F == 1
        win = double(khn>=0);
    elseif F == -1
        ks = (xSourceAct)/norm(xSourceAct);
        win1 = (kh*ks').^2.*((kh*ks')>0);
        v0 = xRef-xSourceAct;
        win = (xSourceAct-x0)*v0';
         win(win<0) = 0;
        %  tukeyTemplate = tukeywin(200,1);
        %  tukeyTemplate = tukeyTemplate(1:end/2);
        % win = interp1(linspace(0,1,length(tukeyTemplate)),tukeyTemplate,win);
        win = win/max(abs(win));
    end

    dref = sqrt(abs( rhoG.*rhoP./(rhoG+F*rhoP )));
    tau = F*rhoP/c;
    amp = 1/(4*pi)*1./rhoP*norm(xSourceAct)*4*pi;
    Dx = -k.*sqrt(F*8*pi./(1i*k)).*(amp.*win.*dref.*1i.*khn.*dl)'.*exp(-1i*omega.*tau');

    % % Give analytical lower freq. limit
    rP0 = min(rhoP);
    rG0 = min(rhoG);
    cG = (1./rG0);
    cP =  (1./rP0);
    K1 = (2*cP.^2+cG.^2)./(F.*cP+cG);
    K2 = +3/4*(F*cP.^3+cG.^3)./(F.*cP+cG).^2;
    omMin = c*max( abs(K1 - K2), 2/Rssd);
    freqAnal(si) =  omMin/2/pi;

    rS = norm(xSourceAct);
    K = 1/4*abs( abs(rS-Rssd)./rS/Rssd- 5./abs(rS-Rssd) );
    omMin2 = c*max( K, pi/Rssd); 
    freqAnal2(si) =  omMin2/2/pi;

    % Estimate synthesized field at xRef and mesure lower freq. limit
    R0 = sqrt(sum((xRef-x0).^2,2));
    Ptransfer = sum( 1/(4*pi)*Dx.*exp(-1i*omega*R0'/c)./R0',2);
    L0 = -1;
    ixs = find(20*log10(abs(Ptransfer))<L0+0.1 & 20*log10(abs(Ptransfer))>L0);
    if ~isempty(ixs)
        freqMeasured(si) = freq(ixs(1));
    else
        freqMeasured(si) = 0;
    end

end
%%
%%
figure
plot(xSource(:,1), freqMeasured)
hold on
plot(xSource(:,1), freqAnal)
%plot(xSource(:,1), freqAnal2)
%semilogy(xSource(:,1), freqAnal2)
%ylim([10,2e3])