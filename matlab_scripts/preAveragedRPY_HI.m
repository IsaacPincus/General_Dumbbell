function Omega = preAveragedRPY_HI(a, H, dQ, sigma)
    ll = max([0 sigma-dQ])+10e-20;
    ul = sigma+dQ-10e-20;
    Q = ll:(ul-ll)/1000:ul;
    
    kB_T = 1.38064852*10^-23*(293.15);
    eta = 9.5*10^-4;
    zeta = 6*pi*eta*a;
    c = (H*dQ^2)/(2*kB_T);
    Jeq = (dQ*(dQ^2+(3+2*c)*sigma^2))/(3+2*c)*beta(1/2,1+c);
    psi_eq = (Q.^2/(Jeq)).*(1-(Q-sigma).^2/dQ^2).^c;
    A = zeros(size(Q));
    B = zeros(size(Q));
    for i=1:length(Q)
        if Q(i)>=2*a
            A(i) = 1+2/3*(a/Q(i))^2;
            B(i) = 1-2*(a/Q(i))^2;
        else
            A(i) = (4/3)*(Q(i)/a) - (3/8)*(Q(i)/a)^2;
            B(i) = (1/8)*(Q(i)/a)^2;
        end
    end
    integrand = (3*a)./(4*zeta*Q).*(A+B/2).*psi_eq;
    
%     figure()
%     plot(Q, integrand);
    
    Omega = trapz(Q, integrand);
end