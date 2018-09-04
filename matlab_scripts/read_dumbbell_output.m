%% Time Data from timestepdata.inp
% Always run this first
dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
delay = dtData(1,4);
size_steps=Nrelax_times/delay+1;

%% Read a single file and plot time-series data at different timesteps

file = 'eta.dat';
rawData = dlmread(file, '');

data = [];
for i=1:length(dt)
    data(:,:,i) = rawData((size_steps+1)*(i-1)+2:(size_steps+1)*i, 1:3);
end

for i=1:size_steps
    step_sorted_data = [dt', reshape(data(i,2:3,:),[2,length(dt)])'];
    [Q(i), dQ(i)] = textra(step_sorted_data, 1, 0.1);
end

figure();
hold on
for i=1:length(dt)
    name = ['dt=',num2str(dt(i))];
    e1 = errorbar(data(:,1,i), data(:,2,i), data(:,3,i),...
        'DisplayName',name,'LineWidth',2);
end
errorbar(data(:,1,i), Q, dQ,...
        'DisplayName','TEXTRA','LineWidth',2);
[h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',16,'Interpreter','latex','Box','off');
hold off

%% Read a list of files and plot data at a specific time
files = ["eta_sr0.1.dat", "eta_sr0.01.dat", "eta_sr0.001.dat"];
xvals = [0.1, 0.01, 0.001]; % must be in the same order as in files
time = 7;

yvals = [];
dyvals = [];
clear Q dQ

for j = 1:length(files)
    rawData = dlmread(files(j), '');

    data = [];
    for i=1:length(dt)
        data(:,:,i) = rawData((size_steps+1)*(i-1)+2:(size_steps+1)*i, 1:3);
    end
    
    t = data(:,1,1);
    [~, timestep] = min(abs(t-time));
    
    step_sorted_data = [dt', reshape(data(timestep,2:3,:),[2,length(dt)])'];
    [Q(j), dQ(j)] = textra(step_sorted_data, 1, 0.25);
end

figure();
hold on
errorbar(xvals, Q, dQ,...
        'DisplayName','TEXTRA','LineWidth',2);
dim = [0.18 0.18 0.7 0.7];
str = ['t = ', num2str(time)];
annotation('textbox',dim,'String',str,'FitBoxToText',...
    'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
hold off

