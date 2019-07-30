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
clearvars -except times chiGsdt

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
Wi = 50;
h = 0.375;
mu1 = 1-h*(1+32/27*h^2);
mu2 = 1-2*h*(1-32/27*h^2);

%N must be even! Order of expansion of spherical harmonics
N = 20;
% tauspan = [0 1];
tauspan = linspace(0,1,150);
A_0 = zeros((N/2+1)*((N/2+1)+1),1);
A_0(1) = 1;

%Perform integration of system of ODEs
[tau, A] = ode15s(@(tau,A) psi_harmonics(tau,A,N,Wi,mu1), tauspan, A_0);

% tau_red = tau(1:ceil(length(tau)/50):end);

%Calculate polymer contribution to viscosity, psi1, psi2
for i = 1:length(tau)
%     time_step = tau_red(i);
%     t = length(tau(tau<time_step))+1;
%     t = tau(i);
    
    ha = A(i,:);
    
    sc(i) = (1/3) - (1/15)*get_A(ha,0,2,0,N) ...
        + (2/5)*get_A(ha,0,2,2,N);
    
    S(i) = 0.5*(3*sc(i) - 1);
    
    eta_p(i) = -(-(6/(5*Wi))*get_A(ha,1,2,2,N) ...
        -(1)/(12*mu2)*((2/5) - (4/35)*get_A(ha,0,2,0,N) + ...
        (2/105)*get_A(ha,0,4,0,N)-16*get_A(ha,0,4,4,N)));
    Psi1_p(i) = -(-(12/(5*Wi^2))*get_A(ha,0,2,2,N)...
        - (8/(Wi*3*mu2))*get_A(ha,1,4,4,N));
    Psi2_p(i) = -((3/Wi^2)*(-(1/5)*get_A(ha,0,2,0,N) ...
        - (2/5)*get_A(ha,0,2,2,N)) ...
        + (1/(Wi*4*mu2))*((8/35)*get_A(ha,1,2,2,N)...
        - (4/7)*get_A(ha,1,4,2,N) -(16/3)*get_A(ha,1,4,4,N)));

    chi_tau(i) = (1/2)*atan(2*eta_p(i)/(Psi1_p(i)*Wi));

    S_tau(i) = 0.5*(3*cos(chi_tau(i))^2-1);
    
    chi_G(i) = (1/2)*atan(get_A(ha,1,2,2,N)/(get_A(ha,0,2,2,N)));
    
    S_G(i) = 0.5*(3*cos(chi_G(i))^2-1);
                    
end

% figure();
% axes1 = gca;
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
% pbaspect([1. 1. 1]);
% xlabel('time ($t/\lambda$)', 'Interpreter', 'latex', 'FontSize', 20')
% ylabel('$S$', 'Interpreter', 'latex', 'FontSize', 20')
% hold on
% e1 = plot(tau, S, 'bo', 'DisplayName','S-parameter','LineWidth',2);
% e1.MarkerFaceColor='b';
% e1.MarkerSize=10;
% e1 = plot(tau, S_tau, 'r>', 'DisplayName','S_tau-parameter','LineWidth',2);
% e1.MarkerFaceColor='r';
% e1.MarkerSize=10;
% e1 = plot(tau, S_G, 'gs', 'DisplayName','S_G-parameter','LineWidth',2);
% e1.MarkerFaceColor='g';
% e1.MarkerSize=10;
% dim = [0.55 0.15 0.3 0.3];
% str = {['$N = $' num2str(N)], ['$\dot{\gamma}\lambda = $' num2str(Wi)]...
%     ,['$h = $' num2str(h)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);
% hold off

figure();
axes1 = gca;
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
pbaspect([1. 1. 1]);
xlabel('time ($t/\lambda$)', 'Interpreter', 'latex', 'FontSize', 20')
ylabel('$S$', 'Interpreter', 'latex', 'FontSize', 20')
hold on
e1 = plot(tau, chi_tau, 'r>', 'DisplayName','chi_tau','LineWidth',2);
e1.MarkerFaceColor='r';
e1.MarkerSize=10;
e1 = plot(tau, chi_G, 'gs', 'DisplayName','chi_g','LineWidth',2);
e1.MarkerFaceColor='g';
e1.MarkerSize=10;
% e1 = plot(tau, sc, 'bh', 'DisplayName','cos^2','LineWidth',2);
% e1.MarkerFaceColor='b';
% e1.MarkerSize=10;
dim = [0.55 0.15 0.3 0.3];
str = {['$N = $' num2str(N)], ['$\dot{\gamma}\lambda = $' num2str(Wi)]...
    ,['$h = $' num2str(h)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);
hold off

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
% 
% hold on
% yyaxis(axes1, 'right')
% axes1.YAxis(2).Color = 'r';
% ylabel('$\frac{\Psi_1}{nkT\lambda^2}$', 'Interpreter', 'latex', 'FontSize', 32')
% e1 = plot(tau,Psi1_p,'rd', 'DisplayName','\Psi_1','LineWidth',2);
% e1.MarkerFaceColor='r';
% e1.MarkerSize=10;
% 
% % yval = 0.5;
% % xval = 0.61;
% % annotation('arrow', [xval xval+0.1], [yval yval], 'LineWidth', 1.5, 'Color', 'k')
% % annotation('arrow', [xval xval], [yval yval-0.1], 'LineWidth', 1.5, 'Color', 'k')
% 
% hold off
% 
% dim = [0.55 0.15 0.3 0.3];
% str = {['$N = $' num2str(N)], ['$\dot{\gamma}\lambda = $' num2str(Wi)]...
%     ,['$h = $' num2str(h)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',16);
% 
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

%% Animation of distribution function
azi_r = -pi:pi/12:pi;
elev_r = -pi/2:pi/12:pi/2;

% azi_m = (azi_r(1:end-1)+azi_r(2:end))/2;
% elev_m = (elev_r(1:end-1)+elev_r(2:end))/2;

[azi_g, elev_g] = meshgrid(azi_r, elev_r);

dist = psi_eq(A(1,:),elev_g,azi_g,N);
[X,Y,Z] = sph2cart(azi_g, elev_g, dist);

k = [0 1 0; 0 0 0; 0 0 0];
minX = -1;
maxX = 1;
minY = -1;
maxY = 1;
minZ = -1;
maxZ = 1;

[Xs, Ys, Zs] = meshgrid(...
    minX:(maxX-minX)/5:maxX,...
    minY:(maxY-minY)/5:maxY,...
    minZ:(maxZ-minZ)/5:maxZ);

for i=1:size(Xs,1)
    for j=1:size(Ys,2)
        for p=1:size(Zs,3)
            fl(:,i,j,p) = k*[Xs(i,j,p);Ys(i,j,p);Zs(i,j,p)];
        end
    end
end

aniFig = figure('units','normalized','outerposition',[0 0 1 1]);
axis equal
hold on
xlabel('x')
ylabel('y')
zlabel('z')
hAni = surf(X,Y,Z);
quiver3(Xs,Ys,Zs,...
    squeeze(fl(1,:,:,:)),squeeze(fl(2,:,:,:)),squeeze(fl(3,:,:,:)),...
    'linewidth', 2, 'MaxHeadSize', 0.5);
view(30,30);
pause(2)
for time=2:size(A,1)
    dist = psi_eq(A(time,:),elev_g,azi_g,N);
    [X,Y,Z] = sph2cart(azi_g, elev_g, dist);
    hAni.XData = X;
    hAni.YData = Y;
    hAni.ZData = Z;
    drawnow limitrate
    pause(0.2)
end
hold off

%% Steady state calculations
clear variables

% FF spring zero-shear viscosity
a_vals = [2.66666666666667E-08 6.66666666666667E-08 1.33333333333333E-07 2E-07 3.33333333333333E-07 5E-07];
for i=1:length(a_vals)
    a = a_vals(i);
    sigma = 1.00E-06;
    dQ = 1E-7;
    H = 8.14996821356E-07;
    kB_T = 1.38064852*10^-23*(293.15);
    eta = 9.5*10^-4;
    c = (H*dQ^2)/(2*kB_T);

    omega = preAveragedRPY_HI(a, H, dQ, sigma);
    zeta = 6*pi*eta*a;

    Q2_eq = (3*dQ^4 + 6*(5+2*c)*dQ^2*sigma^2+(5+2*c)*(3+2*c)*sigma^4)/...
             ((5+2*c)*(dQ^2+(3+2*c)*sigma^2));

    zeta_HI = 1/(1/zeta-omega);

    eta0_p = zeta_HI*Q2_eq/12;

    eta0_p_rodlike(i) = eta0_p/(sigma^2*zeta);
end

% Rodlike viscosity scaling

h_vals = [0 1/8 2/8 3/8];
Wi_vals = logspace(-3,3, 40);
%N must be even!!
N = 40;

for i = 1:length(Wi_vals)
    for j = 1:length(h_vals)
        tic
        h = h_vals(j);
        mu1 = 1-h*(1+32/27*h^2);
        mu2 = 1-2*h*(1-32/27*h^2);
        Wi = Wi_vals(i);
        ha = solve_eq(N, Wi, mu1);
        
        eta(i,j) = -(-(6/(5*Wi))*get_A(ha,1,2,2,N) ...
            -(1)/(12*mu2)*((2/5) - (4/35)*get_A(ha,0,2,0,N) + ...
            (2/105)*get_A(ha,0,4,0,N)-16*get_A(ha,0,4,4,N)));
        Psi1(i,j) = -(-(12/(5*Wi^2))*get_A(ha,0,2,2,N)...
            - (8/(Wi*3*mu2))*get_A(ha,1,4,4,N));
        Psi2(i,j) = -((3/Wi^2)*(-(1/5)*get_A(ha,0,2,0,N) ...
            - (2/5)*get_A(ha,0,2,2,N)) ...
            + (1/(Wi*4*mu2))*((8/35)*get_A(ha,1,2,2,N)...
            - (4/7)*get_A(ha,1,4,2,N) -(16/3)*get_A(ha,1,4,4,N)));
        sc = (1/3) - (1/15)*get_A(ha,0,2,0,N) ...
            + (2/5)*get_A(ha,0,2,2,N);
    
        S(i,j) = 0.5*(3*sc - 1);
        
        chi_tau(i,j) = (1/2)*atan(2*eta(i,j)/Psi1(i,j));
        
        S_tau(i,j) = cos(chi_tau(i,j))^2;
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
e1 = plot(Wi_vals, eta(:,1), 'b--', ...
        'DisplayName','Rod, h = 0','LineWidth',2);
e1.MarkerFaceColor='b';
e1.MarkerSize=14;
e1 = plot(Wi_vals, eta(:,2), 'r--', ...
        'DisplayName','Rod, h = 1/8','LineWidth',2);
e1.MarkerFaceColor='r';
e1.MarkerSize=14;
e1 = plot(Wi_vals, eta(:,3), 'g--', ...
        'DisplayName','Rod, h = 2/8','LineWidth',2);
e1.MarkerFaceColor='g';
e1.MarkerSize=14;
e1 = plot(Wi_vals, eta(:,4), 'm--', ...
        'DisplayName','Rod, h = 3/8','LineWidth',2);
e1.MarkerFaceColor='m';
e1.MarkerSize=14;

FF_strings = ["FF Spring, h=0", "FF Spring, h=1/8", "FF Spring, h=2/8",...
    "FF Spring, h=3/8"];
line_colours = ["b", "r", "g", "m"];

for i=1:length(a_vals)
    line([min(Wi_vals) max(Wi_vals)], [eta0_p_rodlike(i) eta0_p_rodlike(i)],...
        'DisplayName', FF_strings(i), 'color', line_colours(i));
end

ylim([min(eta(:,1))*0.9, max(eta(:,4))*1.1]);

[h,icons,plots,legend_text]=legend({},'Location','southwest','FontSize',16,...
    'Interpreter','latex','Box','off');
xlabel('$\dot{\gamma}\lambda$', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('$\frac{(\eta - \eta_s)}{nkT\lambda}$', 'Interpreter', 'latex', 'FontSize', 32)
% dim = [0.18 0.18 0.7 0.7];
% str = ['t = ', num2str(time)];
% annotation('textbox',dim,'String',str,'FitBoxToText',...
%     'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
hold off

%% S-plot
clear variables

eta = 9.5*10^-4;            %Fluid dynamic viscosity, Pa.s
kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
T = 295.15;                 %Temperature, K
N = 40;

k_vals = linspace(0,3000,10);

for i=1:length(k_vals)
    k = k_vals(i);
    
    L = 800*10^-9;              %Length of dumbbell rod, nm
    a = 200*10^-9;               %Size of beads, nm
    zeta = 6*pi*eta*a;          %Stokes law friction factor, Pa.s.m
    h = (3*a)/(4*L);            %Hydrodynamic interaction parameter
    lambda = (zeta*L^2)/(kB*T);
    Wi = k*lambda;
    mu1 = 1-h*(1+32/27*h^2);
    mu2 = 1-2*h*(1-32/27*h^2);
    
    ha = solve_eq(N, Wi, mu1);
    
    sc = (1/3) - (1/15)*get_A(ha,0,2,0,N) ...
            + (2/5)*get_A(ha,0,2,2,N);
    
    S(i) = 0.5*(3*sc - 1);
    
end

kLD = [780, 1250, 1880, 2500, 3150];
LDr = [0.71, 0.87, 0.96, 1.05, 1.01];

figure();
hold on
plot(k_vals,S, 'go')
plot(kLD, LDr/1.8063, 'rx')
hold off

%% S-plot with lengths
clear variables

eta = 9.5*10^-4;            %Fluid dynamic viscosity, Pa.s
kB = 1.38064852*10^-23;     %Boltzmann constant, m^2 kg s^-2 K^-1
T = 295.15;                 %Temperature, K
N = 40;
k = 3000;
L_vals = linspace(0,1000,10)*10^-9;

for i=1:length(L_vals)
    L = L_vals(i);
    
    L = 800*10^-9;
    a = 10*10^-9;               %Size of beads, nm
    zeta = 6*pi*eta*a;          %Stokes law friction factor, Pa.s.m
    h = (3*a)/(4*L);            %Hydrodynamic interaction parameter
    lambda = (zeta*L^2)/(kB*T);
    Wi = k*lambda;
    mu1 = 1-h*(1+32/27*h^2);
    mu2 = 1-2*h*(1-32/27*h^2);
    
    ha = solve_eq(N, Wi, mu1);
    
    sc = (1/3) - (1/15)*get_A(ha,0,2,0,N) ...
            + (2/5)*get_A(ha,0,2,2,N);
    
    S(i) = 0.5*(3*sc - 1);
    
end

kLD = [780, 1250, 1880, 2500, 3150];
LDr = [0.71, 0.87, 0.96, 1.05, 1.01];

expData = [773.0547550432277, 0.7086614173228347
1257.2046109510088, 0.8645669291338585
1884.005763688761, 0.9842519685039371
2512.9682997118157, 1.0472440944881891
3141.9308357348705, 1.0125984251968505];

kLD = expData(:,1);
LDr = expData(:,2);

figure();
hold on
plot(k_vals,S, 'go')
plot(kLD, LDr/1.8063, 'rx')
hold off

%% Some testing

azi_r = 0:pi/12:2*pi;
elev_r = 0:pi/12:pi;

% azi_m = (azi_r(1:end-1)+azi_r(2:end))/2;
% elev_m = (elev_r(1:end-1)+elev_r(2:end))/2;

[azi_g, elev_g] = meshgrid(azi_r, elev_r);


dist = psi_eq(A(end,:),azi_g,elev_g,N)
% dist = abs(dist/max(max(dist)));

[X,Y,Z] = sph2cart(azi_g, elev_g, dist);
% X = [X;X(1,:)];
% Y = [Y;Y(1,:)];
% Z = [Z;Z(1,:)];

aniFig = figure('units','normalized','outerposition',[0 0 1 1]);
axis equal
hold on
hAni = surf(X,Y,Z);
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
