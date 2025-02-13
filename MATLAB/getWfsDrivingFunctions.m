function [Dx0,dx0,wcMax, omMin] = getWfsDrivingFunctions(model,position, x0, n0, xref, freq, c, AA)
R = sqrt(sum(x0.^2,2));
rhoG = sqrt(sum((x0-xref).^2,2));
dl = mean(R*2*pi/size(x0,1));
omega = 2*pi*freq;
k = omega/c;
switch model
    case 'PW'
        Rs = inf;
        phi0 = atan2(position(2),position(1));
        kh = -[cos(phi0) sin(phi0)];
        khn = n0 * kh';
        win = double(khn>0);
        dref = sqrt(rhoG);
        tau = x0*kh'/c;
        amp = 1;
        tau0 = 0;
        Dx0 = -k.*sqrt(8*pi./(1i*k)).*(amp.*win.*dref.*1i.*khn.*dl)'.*exp(-1i*omega.*(tau' - tau0));
    case 'PS'
        Rs = norm(position);
        rhoP = sqrt(sum( (x0-position).^2,2 ));
        kh = (x0-position)./rhoP;
        khn = sum(n0 .* kh,2);
        F = -1+2*all(khn<=0); % Focused flag
        if F == -1
            win = double(khn>=0);
        elseif F == 1
            ks = (position)/norm(position);
            %win = (kh*ks'+1)/2;
            win = (kh*ks').^2.*((kh*ks')>0);
        end
        dref = sqrt(rhoG.*rhoP./(rhoG-F*rhoP));
        tau = -F*rhoP/c;
        amp = 1/(4*pi)*1./rhoP*norm(position)*4*pi;
        tau0 =  norm(position)/c;
        Dx0 = -k.*sqrt(-F*8*pi./(1i*k)).*(amp.*win.*dref.*1i.*khn.*dl)'.*exp(-1i*omega.*(tau' - tau0));
        %Dx0 = sqrt(-F*8*pi./(1i*k)).*(amp.*win.*dref.*dl.*khn)'.*(-1./rhoP'+1i*k).*exp(-1i*omega.*(tau' - tau0));
        rP0 = min(rhoP);
        rG0 = min(rhoG);

        cG0 = 1./rG0;
        cP0 = 1./rP0;
        K1 = (2*cP0.^2+cG0.^2)./(F.*cP0+cG0);
        K2 = 3/4*(F*cP0.^3+cG0.^3)./(F.*cP0+cG0).^2;

        omMin = c.*abs(K1-K2);
end
Dx0(isnan(Dx0)) = 0;
Dx0(isinf(Dx0)) = 0;

dx1 = (circshift(x0,-1)-x0);
x1 = x0+dx1/2;
dx2 = (circshift(x0,1)-x0);
x2 = x0+dx2/2;
kh1 = (x1-position)./sqrt(sum( (x1-position).^2,2)) ;
kh2 = (x2-position)./sqrt(sum( (x2-position).^2,2));

if AA
    R0 = mean(sqrt(sum(x0.^2,2)));
    [~,ix] = sort(x0*position','descend');
    Nbut = 4;
    dx_dir = max( [sum(kh1.*dx1,2),sum(kh2.*dx2,2)] ,[] ,2);
    ixs_to_mod = ix(1:2);
    G = (x0(ixs_to_mod,:)')\position'/Rs*R0;
    ix0 = find(abs(G-1)<1e-9);
    if ~isempty(ix0)
        ixs_to_mod = ixs_to_mod(ix0);
        G = G(ix0);
    end
    wc = pi.*c./abs(dx_dir);
    wcMax = max(wc.*win);
    [Wc,W] = meshgrid(wc,omega);
    antiAliastingFilters = 1./sqrt( 1 + ( W./Wc ).^(2*Nbut));
    H0 = 1./sqrt(1+(1i*W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ));
    H1 = 1./sqrt(1+(W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ).^(2*Nbut));
    %    G1 = H1 + (1-H1).*sqrt(G').*exp(1i*W(:,ixs_to_mod).*(tau(ixs_to_mod)'-mean(tau(ixs_to_mod))));

    G1 = H1 + (1-H1).*sqrt(G').*exp(1i*W(:,ixs_to_mod).*(tau(ixs_to_mod)'-mean(tau(ixs_to_mod))));
    antiAliastingFilters(:,ixs_to_mod) = H0.*G1;
    equalizationFilters = antiAliastingFilters;
    Dx0 = equalizationFilters .*Dx0;
else
    wcMax = inf;
end
Dx0(isnan(Dx0)) = 0;
Dx0(isinf(Dx0)) = 0;

dx0 = fftshift(ifft(Dx0,(length(freq)-1)*2,1,'symmetric'),1);
end