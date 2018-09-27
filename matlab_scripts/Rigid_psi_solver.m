%This Matlab script is intended to solve for the orientational distribution
%function of a prolate spheriod or rigid dumbbell in shear flow,
%as per McLachlan et al.'s 2013 paper (Calculations of flow-induced
%orientation distributions for analysis of linear dichroism spectroscopy).
%It then spits out the specific S value for this orientational distribution
%function by integrating over the distribution, as well as viscosity and
%first and second normal stress differences.

%This script uses the same non-dimensionalisation (in other words, the same
%definition of lambda) as the rodlike time constant in Larson et al. 2006.
%It also considers hydrodynamic interaction from the RPY tensor, rather
%than the Oseen-Burgers tensor.

% Toggle this to suppress the singular matrix warnings, since they're
% pretty annoying and the results don't appear innaccurate. Be aware that
% it seems you can't set N too high (~>50 or higher). 
warning('off', 'MATLAB:nearlySingularMatrix')

%% Transient behaviour
clear variables

eta = 9.5*10^-4;            %Fluid dynamic viscosity, Pa.s
kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
T = 295.15;                 %Temperature, K
k = 2000;                    %shear rate, s^-1

%here we calculate the dimensionless shear rate, the Weissenberg number Wi
%for prolate spheriods
% a = 400*10^-9;              %semi-diameter parallel to axis of rev, m
% b = 4*10^-9;                %semi-diameter perpendicular to axis of rev, m
% r = a/b;
% Ft = sqrt(r^2-1)/(r^(1/3)*log(r+sqrt(r^2-1)));
% Fr = (4*(r^4-1))/(3*r^2*((2*(2*r^2-1))/(r^(4/3)*Ft))-2);
% lambda = (4*Fr*pi()*eta*a*b^2)/(3*kB*T)*(1/12);
% Wi = k*lambda
% h = 0;
% mu1 = 1;
% mu2 = 1;

% For a bead-rod dumbbell with HI (Stewart and Sorensen 1972, Soc. Rheol):
L = 400*10^-9;              %Length of dumbbell rod, nm
a = 200*10^-9;               %Size of beads, nm
zeta = 6*pi*eta*a;          %Stokes law friction factor, Pa.s.m
h = (3*a)/(4*L);            %Hydrodynamic interaction parameter
lambda = (zeta*L^2)/(kB*T);
Wi = k*lambda;
Wi = 40;
h = 3/8;
mu1 = 1-h*(1+32/27*h^2);
mu2 = 1-2*h*(1-32/27*h^2);

%N must be even! Order of expansion of spherical harmonics
N = 40;
% tauspan = [0 1];
tauspan = linspace(0,1,50);
A_0 = zeros((N/2+1)*((N/2+1)+1),1);
A_0(1) = 1;

%Perform integration of system of ODEs
[tau, A] = ode15s(@(tau,A) psi_harmonics(tau,A,N,Wi,mu1), tauspan, A_0);

% tau_red = tau(1:ceil(length(tau)/50):end);
% Calculate S-parameter
for i = 1:length(tau)
%     time_step = tau(i);
%     t = length(tau(tau<time_step))+1;
%     t = tau(i);
    ha = A(i,:);
    
    sc = (1/3) - (1/15)*get_A(ha,0,2,0,N) ...
        + (2/5)*get_A(ha,0,2,2,N);
    
    S(i) = 0.5*(3*sc - 1);
end

figure();
axes1 = gca;
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
pbaspect([1. 1. 1]);
xlabel('time ($t/\lambda$)', 'Interpreter', 'latex', 'FontSize', 20')
ylabel('$S$', 'Interpreter', 'latex', 'FontSize', 20')
hold on
e1 = plot(tau, S, 'bo', 'DisplayName','S-parameter','LineWidth',2);
e1.MarkerFaceColor='b';
e1.MarkerSize=10;
dim = [0.55 0.15 0.3 0.3];
str = {['$N = $' num2str(N)], ['$\dot{\gamma}\lambda = $' num2str(Wi)]...
    ,['$h = $' num2str(h)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);
hold off

%Calculate polymer contribution to viscosity, psi1, psi2
for i = 1:length(tau)
%     time_step = tau_red(i);
%     t = length(tau(tau<time_step))+1;
%     t = tau(i);
    
    ha = A(i,:);
    eta_p(i) = -(-(6/(5*Wi))*get_A(ha,1,2,2,N) ...
        -(1)/(12*mu2)*((2/5) - (4/35)*get_A(ha,0,2,0,N) + ...
        (2/105)*get_A(ha,0,4,0,N)-16*get_A(ha,0,4,4,N)));
    Psi1_p(i) = -(-(12/(5*Wi^2))*get_A(ha,0,2,2,N)...
        - (4/(Wi*3*mu2))*get_A(ha,1,4,4,N));
    Psi2_p(i) = -((3/Wi^2)*(-(1/5)*get_A(ha,0,2,0,N) ...
        - (2/5)*get_A(ha,0,2,2,N)) ...
        + (1/(Wi*4*mu2))*((8/35)*get_A(ha,1,2,2,N)...
        - (4/7)*get_A(ha,1,4,2,N) -(16/3)*get_A(ha,1,4,4,N)));
end

figure();
set(gcf, 'Position', [500 400 1000 700])
axes1 = gca;
% axes1.XScale = 'log';
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
pbaspect([1. 1. 1]);
ylim([0 max(eta_p)*1.1]);
% hold on
% yyaxis(axes1, 'left')
% axes1.YAxis(1).Color = 'b';
xlabel('time ($t/\lambda$)', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('$\frac{(\eta - \eta_s)}{nkT\lambda}$', 'Interpreter', 'latex', 'FontSize', 32)
e1 = plot(tau,eta_p,'bo', 'DisplayName','Viscosity','LineWidth',2);
e1.MarkerFaceColor='b';
e1.MarkerSize=10;
hold off

% hold on
% yyaxis(axes1, 'right')
% axes1.YAxis(2).Color = 'r';
% ylabel('$\frac{\Psi_1}{nkT\lambda^2}$', 'Interpreter', 'latex', 'FontSize', 32')
% e1 = plot(tau,Psi1_p,'rd', 'DisplayName','\Psi_1','LineWidth',2);
% e1.MarkerFaceColor='r';
% e1.MarkerSize=10;

% yval = 0.5;
% xval = 0.61;
% annotation('arrow', [xval xval+0.1], [yval yval], 'LineWidth', 1.5, 'Color', 'k')
% annotation('arrow', [xval xval], [yval yval-0.1], 'LineWidth', 1.5, 'Color', 'k')

% hold off

dim = [0.55 0.15 0.3 0.3];
str = {['$N = $' num2str(N)], ['$\dot{\gamma}\lambda = $' num2str(Wi)]...
    ,['$h = $' num2str(h)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);

% axes2 = axes('Position',[.4 .25 .37 .37]);
% hold(axes2,'on');
% box(axes2,'on');
% set(axes2,'FontSize',14,'LineWidth',1.5,'TickLength',[0.015 0.025]);
% pbaspect([1. 1. 1]);
% hold on
% xlabel('time ($t/\lambda$)', 'Interpreter', 'latex', 'FontSize', 16)
% ylabel('$\frac{\Psi_2}{nkT\lambda^2}$', 'Interpreter', 'latex', 'FontSize', 30)
% e1 = plot(tau,Psi2_p,'gs', 'DisplayName','\Psi_2','LineWidth',2);
% e1.MarkerFaceColor='g';
% e1.MarkerSize=5;
% hold off

%% Steady state calculations
clear variables
h_vals = [0 1/8 2/8 3/8];
Wi_vals = logspace(-1,3, 20);
%N must be even!!
N = 40;

for i = 1:length(Wi_vals)
    for j = 1:length(h_vals)
        tic
        h = h_vals(j);
        mu1 = 1-h*(1+32/27*h^2);
        mu2 = 1-h*(1-32/27*h^2);
        Wi = Wi_vals(i);
        ha = solve_eq(N, Wi, mu1);
        
        eta(i,j) = -(-(6/(5*Wi))*get_A(ha,1,2,2,N) ...
            -(1)/(12*mu2)*((2/5) - (4/35)*get_A(ha,0,2,0,N) + ...
            (2/105)*get_A(ha,0,4,0,N)-16*get_A(ha,0,4,4,N)));
        Psi1(i,j) = -(-(12/(5*Wi^2))*get_A(ha,0,2,2,N)...
            - (4/(Wi*3*mu2))*get_A(ha,1,4,4,N));
        Psi2(i,j) = -((3/Wi^2)*(-(1/5)*get_A(ha,0,2,0,N) ...
            - (2/5)*get_A(ha,0,2,2,N)) ...
            + (1/(Wi*8*mu2))*((8/35)*get_A(ha,1,2,2,N)...
            - (4/7)*get_A(ha,1,4,2,N) -(16/3)*get_A(ha,1,4,4,N)));
        sc = (1/3) - (1/15)*get_A(ha,0,2,0,N) ...
            + (2/5)*get_A(ha,0,2,2,N);
    
        S(i,j) = 0.5*(3*sc - 1);
    end
end

figure();
hold on
axes1 = gca;
fsize=20;
pbaspect([1. 1. 1]);
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
% Comment out the below 2 lines
% as needed.
axes1.XScale='log';
axes1.YScale='log';
e1 = plot(Wi_vals, eta(:,1), 'bd-', ...
        'DisplayName','h = 0','LineWidth',2);
e1.MarkerFaceColor='b';
e1.MarkerSize=14;
e1 = plot(Wi_vals, eta(:,end), 'ro-', ...
        'DisplayName','h = 3/8','LineWidth',2);
e1.MarkerFaceColor='r';
e1.MarkerSize=14;
[h,icons,plots,legend_text]=legend({},'Location','southwest','FontSize',16,...
    'Interpreter','latex','Box','off');
xlabel('$\dot{\gamma}\lambda$', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('$\frac{(\eta - \eta_s)}{nkT\lambda}$', 'Interpreter', 'latex', 'FontSize', 32)
% dim = [0.18 0.18 0.7 0.7];
% str = ['t = ', num2str(time)];
% annotation('textbox',dim,'String',str,'FitBoxToText',...
%     'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
hold off

%% Functions

function A = get_A(Avector, i, n, m, N)
    %Since my internal numbering for the A vector is very convoluted, I've
    %written this function to get A(i,n,m) from the Avector.
    %Numbering is A(1) = A000, A(2) = A010, A(3) = A002, A(4) = A012 etc...
    %where it's Anim (confusingly...). So even indices is i=1, odd indices
    %is i=0, then you count over m first then n in the order of the sum.
    
    if i==1
        j = 2;
        for q = 0:2:N
            for p= 0:2:q
                if (p==m)&&(q==n)
                    A = Avector(j);
                    return;
                else
                    j = j + 2;
                end
            end
        end
    else
        j = 1;
        for q = 0:2:N
            for p= 0:2:q
                if (p==m)&&(q==n)
                    A = Avector(j);
                    return;
                else
                    j = j + 2;
                end
            end
        end
    end
end

function psi = psi_eq(harmonics, theta, phi, N)
    %This function calculates the value of the orientational distribution
    %function over theta and phi, given the harmonics
    %and the expansion order of the harmonics
    
    i = 1;
    psi = zeros(size(theta));
    for n = 0:2:N
        for m = 0:2:n
            Pmn = legendre(n,cos(theta));
            if n~=0
                Pmn = permute(Pmn, [2 3 1]);
            end
            psi = psi + Pmn(:,:,m+1).*...
                (harmonics(i)*cos(m*phi) + harmonics(i+1)*sin(m*phi));
            i = i+2;
        end
    end
    psi = psi/(4*pi());
    
end

function A = solve_eq(N, Wi, mu)
    %This solves the same system of linear equations as the
    %psi_harmonics function, only at steady state (and far faster)
    %N must be even!!!
    %Numbering is A(1) = A000, A(2) = A010, A(3) = A002, A(4) = A012 etc...
    
    size = (N/2+1)*((N/2+1)+1);
    Am = zeros(size,size);
    
    i = 1;
    for q = 0:2:N
        for p = 0:2:q
            j = 2;
            for n = 0:2:N
                for m= 0:2:n
                    Am(i,j) = -Wi*amnpq(m,n,p,q);
                    j = j + 2;
                end
            end
            Am(i,i) = -(q*(q+1))*mu*2;
            i = i+1;
            
            j = 1;
            for n = 0:2:N
                for m= 0:2:n
                    if p == 0
                        Am(i,j) = 0;
                        continue;
                    end
                    Am(i,j) = +Wi*amnpq(m,n,p,q);
                    j = j + 2;
                end
            end
            Am(i,i) = -(q*(q+1))*mu*2;
            if p == 0
                Am(i,i) = 1;
            end
            i = i+1;
        end
    end
    Am(1,1) = 1;
    x = zeros(size,1);
    x(1) = 1;
    A = Am\x;
end

function Adash = psi_harmonics(t, y, N, Wi, mu)
    %This function will essentially return the system of coupled ODEs given
    %in equations B.2a and B.2b in McLachlan et al. This should be able to 
    %account for an arbitrary N, the paper used N=12
    %N must be even!!!
    %Numbering is y(1) = A000, y(2) = A010, y(3) = A002, y(4) = A012 etc...
    
    i = 1;
    Adash = zeros((N/2+1)*((N/2+1)+1),1);
    
    for q = 0:2:N
        for p = 0:2:q
            j = 2;
            aA_sum = 0;
            for n = 0:2:N
                for m= 0:2:n
                    aA_sum = aA_sum + y(j)*amnpq(m,n,p,q);
                    j = j + 2;
                end
            end
            Adash(i) = -(q*(q+1))*2*mu*y(i) - Wi*aA_sum;
            i = i+1;
            
            j = 1;
            aA_sum = 0;
            for n = 0:2:N
                for m= 0:2:n
                    aA_sum = aA_sum + y(j)*amnpq(m,n,p,q);
                    j = j + 2;
                end
            end
            Adash(i) = -(q*(q+1))*2*mu*y(i) + Wi*aA_sum;
            if p == 0
                Adash(i) = 0;
            end
            i = i+1;
            
        end
    end
    Adash(1) = 0;
end

function a = amnpq(m,n,p,q)
    %This function implements the a^{mp}_{nq} in table B.1 in McLachlan,
    %used in equations B.2 in the same paper
    
    if m==0
        kr = 1;
    else
        kr = 0;
    end
    
    if     (p == m-2)&&(q == n-2)
        a = ((n-2)*factorial(n+m)*(1-kr))/...
              (4*(2*n+1)*(2*n-1)*factorial(n+m-4));
          
    elseif (p == m-2)&&(q == n)
        a = (3*factorial(n-m+2)*factorial(n+m)*(1-kr))/...
              (4*(2*n-1)*(2*n+3)*factorial(n+m-2)*factorial(n-m));
          
    elseif (p == m-2)&&(q == n+2)
        a = -((n+3)*factorial(n-m+4)*(1-kr))/...
              (4*(2*n+1)*(2*n+3)*factorial(n-m));
          
    elseif (p == m)&&(q == n)
        a = -m/2;
        
    elseif (p == m+2)&&(q == n-2)
        a = -((n-2)*(1+kr))/(4*(2*n+1)*(2*n-1));
        
    elseif (p == m+2)&&(q == n)
        a = -(3*(1+kr))/(4*(2*n+3)*(2*n-1));
        
    elseif (p == m+2)&&(q == n+2)
        a = ((n+3)*(1+kr))/(4*(2*n+3)*(2*n+1));
        
    else 
        a = 0;
        
    end
end

% Some extra junk in case you don't trust my derivations and want to check
% the integrals explicitly:
% 
% psi_S_int = @(theta,phi,time)...
%     sin(theta).^3.*cos(phi).^2.*psi_t(tau, A, time, theta, phi, N);
% 
% psi_eta_int = @(theta,phi,time) ...
%     -(-(3/(2*Wi))*(1-h)/(1-2*h)*sin(theta).^2.*sin(2*phi) ...
%      -(3/2)*sin(theta).^4.*sin(2*phi).^2) ...
%      .*sin(theta).*psi_t(tau, A, time, theta, phi, N);
%  
% psi_Psi1_int = @(theta,phi,time) ...
%     -(-(3/Wi^2)*(1-h)/(1-2*h)*sin(theta).^2.*cos(2*phi) ...
%      -(3/Wi)*sin(theta).^4.*sin(2*phi).*cos(2*phi)) ...
%      .*sin(theta).*psi_t(tau, A, time, theta, phi, N);
% 
% psi_Psi2_int = @(theta,phi,time) ...
%     -(-(3/Wi^2)*(1-h)/(1-2*h)*(sin(theta).^2.*sin(phi).^2 -cos(theta).^2)...
%      -(3/Wi)*(sin(theta).^2.*sin(phi).^2 -cos(theta).^2)...
%      .*sin(theta).^2.*sin(2*phi)) ...
%      .*sin(theta).*psi_t(tau, A, time, theta, phi, N);
%  
% psi_int = @(theta,phi,time)...
%     sin(theta).*psi_t(tau, A, time, theta, phi, N);

% eta_p(i) = integral2(@(theta, phi) psi_eta_int(theta, phi, time_step),...
%     0, pi(), 0, 2*pi());
% Psi1(i) = integral2(@(theta, phi) psi_Psi1_int(theta, phi, time_step),...
%     0, pi(), 0, 2*pi());
% Psi2(i) = integral2(@(theta, phi) psi_Psi2_int(theta, phi, time_step),...
%     0, pi(), 0, 2*pi());
% sc = integral2(@(theta, phi) psi_S_int(theta, phi, time_step),...
%     0, pi(), 0, 2*pi());
% 
% S(i) = 0.5*(3*sc - 1);

% function psi = psi_t(time_vector, harmonics, time, theta, phi, N)
%     %This function calculates the value of the orientational distribution
%     %function at a particular time, theta and phi, given the harmonics
%     %and the expansion order of the harmonics
%     
%     %I should technically do a linear interpolation...
%     t = length(time_vector(time_vector<time))+1;
%     
%     psi = psi_eq(harmonics(t,:), theta, phi, N);
% end
