% a = -20010.0000000000
% b = 100199999.899999
% c= -999999998.990000
% 
% a^2
% b
% a^3
% a*b
% c
% 
% disp(['a^2-3b=', num2str(a^2-3*b)])

% Q = 0.85:0.0001:1.15;
% p = (Q-1.10008325928356).*(Q-0.899916740600020).*(Q-1000000.00000000);
% 
% plot(Q, p)
% 
% cubic = 1.09999999216552;
% newton= 1.09999999749973;

data = dlmread('failures.dat');

% dt10dat = data(data(:,4)==10,:);
% size(dt10dat);
% dt1dat = data(data(:,4)==1,:);
% dt01dat = data(data(:,4)==0.1,:);
% dt001dat = data(data(:,4)==0.01,:);
% dt0001dat = data(data(:,4)==0.001,:);
% 
% a_e6_dat = data(data(:,3)==1e-6,:);
% a_e5_dat = data(data(:,3)==1e-5,:);
% a_e4_dat = data(data(:,3)==1e-4,:);
% a_e3_dat = data(data(:,3)==1e-3,:);
% a_e2_dat = data(data(:,3)==1e-2,:);
% a_e1_dat = data(data(:,3)==1e-1,:);
% a_1_dat = data(data(:,3)==1e-0,:);
% a_10_dat = data(data(:,3)==1e1,:);
% a_100_dat = data(data(:,3)==1e2,:);

m = 3;
for i=-3:2
    size(data(data(:,m)<(10^(i)+1e-8)&data(:,m)>(10^(i)-1e-8),:))
end

% Q0vals = data(:,1);
% alphavals = data(:,3);

% Q0vals./alphavals

Q0 = 1000;
dt = 0.01;
alpha = logspace(-6,2,9);
Y = logspace(-4,8,13);

data = data(data(:,4)==dt, :);
data = data(data(:,1)==Q0, :);
data = data(:,2:3);

figure();
plot(data(:,1), data(:,2), 'rx');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([min(Y), max(Y)]);
ylim([min(alpha), max(alpha)]);
xlabel('Y')
ylabel('\alpha')
title_str = ['Q0=', num2str(Q0), ' dt=', num2str(dt)];
title(title_str)

% figure();
% loglog(
% 
% figure
% x = [0 0 0 0 1 0 1 1 0 0 1 0 0 0];
% x(x==0) = NaN; % Changing the value of x. If you need it preserved, save to a new variable.
% h = plot(1:numel(x),x,'k.');
% set(h,'MarkerSize',24)
% xlim([0.5 numel(x)+0.5])



