clear
close all
Lx = 50;
dx = 0.1;

x0 = (-Lx:dx:Lx);
y0 = 0;
fu = 2.5e3;
Nf = 5e3;
freq = logspace(-2,log10(fu),Nf)';

%freq = 1e3;
c = 343;
omega = 2*pi*freq;
k = omega/c;

dx2 = 0.1;
xs = (-0:5*dx2:0)';
ys = (-2:dx2:2)';
[xSource,ySource] = meshgrid(xs,ys);
xSource = [xSource(:),ySource(:)];
%%
for si = 1 : size(xSource)
    si
    xSourceAct = xSource(si,:);
    xRef = [0,2];

    F = -sign(xSourceAct(2));

    rP = sqrt( (x0-xSourceAct(1)).^2 + (y0-(xSourceAct(2))).^2 );
    rG = sqrt( (x0-xRef(1)).^2 + (y0-xRef(2)).^2 );

    % Linear reference
    dref = abs(xRef(2)./(xRef(2)-xSourceAct(2)));
    %dref = abs(1./abs(F./rP+1./rG));
    Dx = -1/4/pi*sqrt(8*pi/1i./k).*sqrt(dref).*1i.*k.*xSourceAct(2)./sqrt(rP).*exp(-F*1i*k.*rP)./rP;

    P = zeros(length(freq),1);
    R0 = sqrt( (x0-xRef(1)).^2 + (y0-xRef(2)).^2 );
    for n = 1 : size(Dx,2)
        P = P + 1/4/pi*dx*exp(-1i*k*R0(n))./R0(n).*Dx(:,n)*norm(xRef-xSourceAct)*4*pi;
    end

    xi = line_intersection(xSourceAct', xRef', [x0(1),y0]', [x0(end),y0]')';
    [~, stat_ix] = min(abs(xi(1)-x0));

    rP0 = sqrt( (xi(1)-xSourceAct(1)).^2 + (xi(2)-(xSourceAct(2))).^2 );
    rG0 = sqrt( (xi(1)-xRef(1)).^2 + (xi(2)-xRef(2)).^2 );
    cG0 = 1./rG0;
    cP0 = 1./rP0;
    K1 = (2*cP0.^2+cG0.^2)./(F.*cP0+cG0);
    K2 = 3/4*(F*cP0.^3+cG0.^3)./(F.*cP0+cG0).^2;

    omMin = c.*abs(K1-K2);
    freqAnal(si) = omMin / 2 / pi;

    L0 = -1;
    ixs = find(20*log10(abs(P))<L0+0.1 & 20*log10(abs(P))>L0-0.1);
    if ~isempty(ixs)
        freqMeasured(si) = freq(ixs(1));
    else
        freqMeasured(si) = 0;
    end

end
%%
freqMapMeasured = reshape(freqMeasured,size(ys,1),size(xs,1));
freqMapAnal = reshape(freqAnal,size(ys,1),size(xs,1));
% figure
% subplot(1,3,1)
% freqMapMeasured(ys == 0,:) = nan;
% surf(xs,ys,(freqMapMeasured));
% xlabel('x_s [m]')
% ylabel('y_s [m]')
% title('Measured low frequency limit')
% zlim([0,600])
% 
% subplot(1,3,2)
% freqMapAnal(ys == 0,:) = nan;
% surf(xs,ys,(freqMapAnal));
% xlabel('x_s [m]')
% ylabel('y_s [m]')
% title(['Analyitical low frequency limit'])
% zlim([0,600])
% 
% subplot(1,3,3)
plot(ys,freqMapMeasured(:,ceil(end/2)))
hold on
plot(ys,freqMapAnal(:,ceil(end/2)))
xlabel('x_s [m]')
ylabel('f_{LF} [Hz]')
grid on
title(['Low frequency limit over x_s = 0'])
legend('f_{measured}',['f_{analytical} '])
