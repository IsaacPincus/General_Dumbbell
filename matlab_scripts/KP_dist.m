% This script should calculate the distribution function of KP chain close
% to the rigid-rod limit, using the expressions from 'Radial Distribution
% Function of Semiflexible Polymers', Wilhelm and Frey, Phys. Rev. Let.
% 1996. The key equations are (2) and (3) in the referenced paper. These
% are in units of energy as kB*T and lengths as L. 

kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
T = 295.15;                 %Temperature, K
kT = kB*T;

kappa = 2;
L = 800*10^-9;
lp = 1200*10^-9/L;
% t = 80;
% lp = 1/t;               % persistence length
% t = 1/lp;               % L/lp, so 1/lp in dimensionless form
r = 0:0.001:2;          % Rod lengths
G = zeros(size(r));     % Distribution function
N = 10000;                 % Series in 

for i=1:length(r)
    if r(i) > 1-0.001
        G(i) = 0;
    elseif lp*(1-r(i)) >= 0.2
        for k=1:N
            G(i) = G(i) + pi^2*k^2*(-1)^(k+1)*exp(-lp*pi^2*k^2*(1-r(i)));
        end
        G(i) = lp/(2*pi)*G(i);        
    else 
        for l=1:N
            G(i) = G(i) + 1/(lp*(1-r(i)))^(3/2)*...
                          exp(-(l-0.5)^2/(lp*(1-r(i))))*...
                          (4*((l-0.5)/sqrt(lp*(1-r(i))))^2-2);
        end
        G(i) = lp/(8*pi^(3/2))*G(i);
    end
end

% Normalisation
% G = G/trapz(r,G);

% Changing to true units
r = r*L;
% G = G/(800*10^-9);

% FF chain distribution function
% Hookean Units
% Jeq = (1/(alpha+3)+sigma.^2/alpha)*beta(1/2,(alpha+2)/2)*alpha^1.5;
% psiQ = Q.^2.*(1-(Q-sigma).^2/alpha).^(alpha/2)./Jeq;

% Rodlike Units
% H = 2;
% s = 0.6;
% Q = 1-s+0.00001:0.001:1+s-0.00001;
% Jeq = (1/(3+H/s^2)+1/s^2)*beta(1/2,1+H/(2*s^2))*s^3;
% psiQ = (Q).^2.*(1-((Q)-1).^2/s^2).^(H/(2*s^2))/Jeq;
% 
% % Full units
% dQ = 2000*10^-9;
% sigma = 400*10^-9;
% H = 2*kT/sigma^2;
% Q = sigma-dQ:sigma/1000:sigma+dQ;
% h = H*dQ^2/(2*kT);
% Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
% psiQ = Q.^2.*(1-(Q-sigma).^2/dQ^2).^h/Jeq;

% figure();
% hold on
% plot(r, (4*pi*r.^2.*G)/trapz(r, (4*pi*r.^2.*G)));
% plot(Q, psiQ);
% hold off

% inputs0(1) = H;
% inputs0(2) = dQ;
% inputs0(3) = sigma;
% options = optimset('MaxFunEvals', 10^4);
% outputs = fminsearch(@(inputs) diff_G_psi(inputs, 4*pi*r.^2.*G, r),...
%     inputs0, options);
% dQ = outputs(2);
% sigma = outputs(3);
% H = outputs(1);
% Q = sigma-dQ:sigma/1000:sigma+dQ;
% h = H*dQ^2/(2*kT);
% Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
% psiQ = Q.^2.*(1-(Q-sigma).^2/dQ^2).^h/Jeq;

% figure();
% hold on
% plot(r, (4*pi*r.^2.*G)/trapz(r, (4*pi*r.^2.*G)));
% plot(Q, psiQ);
% hold off

% clear outputs inputs0
% dQ = 100E-9;
% % sigma = 1E-10;
% % H = 2.86125015519553E-06;
% H = 1.5E-8;
% inputs0(1) = H;
% inputs0(2) = dQ;
% options = optimset('MaxFunEvals', 10^4);
% outputs = fminsearch(@(inputs) diff_G_psi_const(inputs, 4*pi*r.^2.*G, r,L),...
%     inputs0, options)
% dQ = outputs(2);
% H = outputs(1);
% sigma = L-dQ;
% % dQ = dQ*5;
% Q = sigma-dQ:sigma/1000:sigma+dQ;
% h = H*dQ^2/(2*kT);
% Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
% psiQ = Q.^2.*(1-(Q-sigma).^2/dQ^2).^h/Jeq;

clear outputs inputs0
dQ = 0.1*L;
sigma = 0.9*L;
% H = 2.86125015519553E-06;
H = 5*kT/dQ^2;
options = optimset('MaxFunEvals', 10^4);
outputs = fminsearch(@(inputs) diff_G_psi_dQ_only(inputs, 4*pi*r.^2.*G, r,L),...
    dQ, options)
dQ = outputs;
sigma = L-dQ;
H = 5*kT/dQ^2;
Q = sigma-dQ:sigma/1000:sigma+dQ;
h = H*dQ^2/(2*kT);
Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
psiQ = Q.^2.*(1-(Q-sigma).^2/dQ^2).^h/Jeq;


figure();
hold on
plot(r, (4*pi*r.^2.*G)/trapz(r, (4*pi*r.^2.*G)));
plot(Q, psiQ);
legend('WLC Distribution', 'FENE Distribution', 'Location', 'NW');
xlabel('Length (m)');
ylabel('\psi');
xlim([0 sigma+dQ]);
dim = [0.55 0.15 0.3 0.3];
str = {['$H = $' num2str(H) ' N/m'], ['$dQ = $' num2str(dQ) ' m'],...
    ['$\sigma = $' num2str(sigma) ' m'], ['$L = $' num2str(L) ' m'],...
    ['$l_p = $' num2str(lp*L) ' m']};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);
hold off

function f = diff_G_psi(inputs, G, r)
    kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
    T = 295.15;                 %Temperature, K
    kBT = kB*T;
    H = inputs(1);
    dQ = inputs(2);
    sigma = inputs(3);
    h = H*dQ^2/(2*kBT);
    Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
    
    % Make sure G is normalised
    G = G/trapz(r, G);
    
    psiQ = zeros(size(G));
    
    for i=1:length(r)
        if (r(i)>(sigma-dQ))&&(r(i)<(sigma+dQ))
            psiQ(i) = r(i).^2.*(1-(r(i)-sigma).^2/dQ^2).^h/Jeq;
        else
            psiQ(i) = 0;
        end
    end
    
    f = sum((psiQ-G).^2);
    
end

function f = diff_G_psi_const(inputs, G, r, L)
    kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
    T = 295.15;                 %Temperature, K
    kBT = kB*T;
    H = inputs(1);
    dQ = inputs(2);
    sigma = L-dQ;
%     if sigma<0
%         sigma=0;
%     end
    h = H*dQ^2/(2*kBT);
    Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
    
    % Make sure G is normalised
    G = G/trapz(r, G);
    
    psiQ = zeros(size(G));
    
    for i=1:length(r)
        if (r(i)>(sigma-dQ))&&(r(i)<(sigma+dQ))
            psiQ(i) = r(i).^2.*(1-(r(i)-sigma).^2/dQ^2).^h/Jeq;
        else
            psiQ(i) = 0;
        end
    end
    
    f = sum((psiQ-G).^2);
    
end

function f = diff_G_psi_dQ_only(dQ, G, r, L)
    kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
    T = 295.15;                 %Temperature, K
    kBT = kB*T;
    sigma = L-dQ;
    H = 5*kBT/dQ^2;
%     if sigma<0
%         sigma=0;
%     end
    h = H*dQ^2/(2*kBT);
    Jeq = (dQ^3/(3+2*h)+sigma^2*dQ)*beta(0.5,1+h);
    
    % Make sure G is normalised
    G = G/trapz(r, G);
    
    psiQ = zeros(size(G));
    
    for i=1:length(r)
        if (r(i)>(sigma-dQ))&&(r(i)<(sigma+dQ))
            psiQ(i) = r(i).^2.*(1-(r(i)-sigma).^2/dQ^2).^h/Jeq;
        else
            psiQ(i) = 0;
        end
    end
    
    f = sum((psiQ-G).^2);
    
end

