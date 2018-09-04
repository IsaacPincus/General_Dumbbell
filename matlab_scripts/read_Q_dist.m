% clear variables
Q0 = 10;
b = 1;

width = 2*sqrt(b)/(100-1);
Q = (Q0-sqrt(b)):width:(Q0+sqrt(b));
Q(1) = Q(1) + 0.00000001;
Q(end) = Q(end) - 0.00000001;

Jeq = (1/(b+3)+Q0.^2/b)*beta(1/2,(b+2)/2)*b^1.5;
psiQ = Q.^2.*(1-(Q-Q0).^2/b).^(b/2)./Jeq;
test2 = cum_trapz(Q,psiQ);
test2 = test2/test2(end);

rands = rand(1000,1);
for i=1:length(rands)
    dist(i) = inverse_extrap(Q,test2,rands(i));
end
example = inverse_extrap(Q,test2,0.99);

bins = 30;

edges_azi = -pi:2*pi/bins:pi;

tic
data = ones(3, 10^7);
fileID = fopen('Q_dist_output_noHI.dat', 'r');
dummy = fgets(fileID);
data = fscanf(fileID, '%f %f %f', size(data))';
toc
[azi, elev, Ql] = cart2sph(data(:,1), data(:,2), data(:,3));
[Qldist_noHI, ~] = histcounts(Ql, bins, 'Normalization', 'pdf');
[azidist_noHI, ~] = histcounts(azi, edges_azi, 'Normalization', 'pdf');
[elevdist_noHI, ~] = histcounts(elev, bins, 'Normalization', 'pdf');


tic
data = ones(3, 10^7);
fileID = fopen('Q_dist_output_withHI.dat', 'r');
dummy = fgets(fileID);
data = fscanf(fileID, '%f %f %f', size(data))';
toc
[azi, elev, Ql] = cart2sph(data(:,1), data(:,2), data(:,3));
[Qldist_withHI, edges_Ql] = histcounts(Ql, bins, 'Normalization', 'pdf');
[azidist_withHI, ~] = histcounts(azi, edges_azi, 'Normalization', 'pdf');
[elevdist_withHI, edges_elev] = histcounts(elev, bins, 'Normalization', 'pdf');
middles_Ql = diff(edges_Ql)/2+edges_Ql(1:end-1);
middles_azi = diff(edges_azi)/2+edges_azi(1:end-1);
middles_elev = diff(edges_elev)/2+edges_elev(1:end-1);

figure();
axes1 = gca;
axes1.XScale = 'log';
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',14,'LineWidth',1.5,'TickLength',[0.015 0.025]);
pbaspect([1. 1. 1]);
xlim([Q0-sqrt(b) Q0+sqrt(b)])
hold on
plot(Q, test2, 'k-', 'DisplayName','Analytical CDF','LineWidth',1.4);
plot(Q, psiQ, 'k--', 'DisplayName','Analytical PDF','LineWidth',1.4);
e1 = plot(middles_Ql, Qldist_noHI, 'ro', 'DisplayName','Simulation, $h^*=0$','LineWidth',2);
e1.MarkerFaceColor='r';
e1.MarkerSize=7;
e1 = plot(middles_Ql, Qldist_withHI, 'bd', 'DisplayName','Simulation, $h^*=1.5$','LineWidth',2);
% e1.MarkerFaceColor='r';
e1.MarkerSize=7;
hold off
[h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',16,'Interpreter','latex','Box','off');
xlabel('$Q_{eq}$','FontSize',14,'Interpreter','latex');
ylabel('Probability Density','FontSize',14,'Interpreter','latex');
dim = [0.55 0.15 0.3 0.3];
str = {'$\Delta t = 0.001$', '$\mathcal{O}(10^7)$ trajectories'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',14);

figure();
axes1 = gca;
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',14,'LineWidth',1.5,'TickLength',[0.015 0.025]);
pbaspect([1. 1. 1]);
xlim([-pi()/2 pi()/2]);
ylim([0 0.7]);
hold on
e1 = plot(middles_elev, elevdist_noHI, 'ro', 'DisplayName','Simulation, $h^*=0$','LineWidth',2);
e1.MarkerFaceColor='r';
e1.MarkerSize=5;
e1 = plot(middles_elev, elevdist_withHI, 'bd', 'DisplayName','Simulation, $h^*=1.5$','LineWidth',2);
e1.MarkerSize=5;
plot(middles_elev, cos(middles_elev)/2, 'k--', 'DisplayName','Equilibrium result','LineWidth',1.5);
hold off
[h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',16,'Interpreter','latex','Box','off');
xlabel('$\phi$','FontSize',14,'Interpreter','latex');
ylabel('Probability Density','FontSize',14,'Interpreter','latex');
dim = [0.55 0.15 0.3 0.3];
str = {'$\Delta t = 0.001$', '$\mathcal{O}(10^7)$ trajectories'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',14);

axes2 = axes('Position',[.57 .22 .25 .25]);
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',14,'LineWidth',1.5,'TickLength',[0.015 0.025]);
pbaspect([1. 1. 1]);
xlim([-pi() pi()]);
ylim([0 0.2]);
hold on
e1 = plot(middles_azi, azidist_noHI, 'ro', 'DisplayName','Simulation, $h^*=0$','LineWidth',2);
e1.MarkerFaceColor='r';
e1.MarkerSize=5;
e1 = plot(middles_azi, azidist_withHI, 'bd', 'DisplayName','Simulation, $h^*=1.5$','LineWidth',2);
e1.MarkerSize=5;
plot([-pi() pi()], [1/(2*pi) 1/(2*pi)], 'k--', 'DisplayName','Equilibrium result','LineWidth',1.5);
hold off
% [h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',16,'Interpreter','latex','Box','off');
xlabel('$\theta$','FontSize',14,'Interpreter','latex');
% ylabel('Probability Density','FontSize',14,'Interpreter','latex');
% dim = [0.55 0.15 0.3 0.3];
% str = {'$\Delta t = 0.001$', '$\mathcal{O}(10^7)$ trajectories'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',14);

function f = cum_trapz(x, fx)
    temp = zeros(size(x));
    for k = 2:length(x)
        temp(k) = temp(k-1) + (fx(k) + fx(k-1))*(x(k) - x(k-1))/2;
    end
    f = temp;
end

function f = inverse_extrap(x,fx,fxval)
    for i=1:length(x)
        if (fx(i) > fxval)
            f = (fxval-fx(i-1))/(fx(i)-fx(i-1))*(x(i)-x(i-1)) + x(i-1);
            break
        end
    end
end