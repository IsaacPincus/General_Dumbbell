%% Time Data from timestepdata.inp
% Always run this first
dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
delay = dtData(1,4);
size_steps=Nrelax_times/delay+1;
inpdata = dlmread('inputparameters.inp', '');
sigma = inpdata(4);
alpha = inpdata(2);
%set Rodlike_opt to true for Rodlike units, false for Hookean
Rodlike_opt = true;

%% Read a single file and plot time-series data at different timesteps

file = 'eta_sr0.0004.dat';
rawData = dlmread(file, '');

data = nan(size_steps,3,length(dt));
pos = 0;
for i=1:length(dt)
    if dt(i)>=delay
        size_steps=Nrelax_times/dt(i);
    else
        size_steps=Nrelax_times/delay+1;
    end
    pos = pos + size_steps + 1;
    data(1:size_steps,1:3,i) = rawData(pos-(size_steps+1)+2:pos, 1:3);
end

for i=1:size_steps
    step_sorted_data = [dt', reshape(data(i,2:3,:),[2,length(dt)])'];
    [Q(i), dQ(i)] = textra(step_sorted_data, 0, 0.1);
end

if(Rodlike_opt)
    data(:,1,:) = data(:,1,:)/(4*sigma^2);
    %for Viscosity data
    data(:,2:3,:) = data(:,2:3,:)/(4*sigma^2);
    Q = Q/(4*sigma^2);
    dQ = dQ/(4*sigma^2);
    %for Psi1, Psi2 data
%     data(:,2:3,:) = data(:,2:3,:)/(16*sigma^4);
%     Q = Q/(16*sigma^4);
%     dQ = dQ/(16*sigma^4);
    
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
xvals = [0.0004, 0.0012, 0.004, 0.012,...
    0.04, 0.12, 0.4, 1.2, 4, 12, 40];
% xvals = dlmread('srvals.dat', '');
% tmaxvals = dlmread('tvals.dat','');

% If you've used a different naming convention, you can just set the 
% files variable as an array of strings directly
prefix = "eta_sr";
suffix = ".dat";
for i=1:length(xvals)
    str = sprintf('%.15f ',xvals(i));
    if floor(xvals(i))==xvals(i)
        str = regexprep(str, '\.[0]+ ', '');
    else
        str = regexprep(str, '[0]+ ', '');
    end
    files(i) = strcat(prefix, str, suffix);
end

% time = 7;

yvals = [];
dyvals = [];
clear Q dQ

for j = 1:length(files)
    rawData = dlmread(files(j), '');
    
%     Nrelax_times = tmaxvals(j);
    data = nan(size_steps,3,length(dt));
    pos = 0;
    for i=1:length(dt)
        if dt(i)>delay
            size_steps=Nrelax_times/dt(i);
        else
            size_steps=Nrelax_times/delay+1;
        end
        pos = pos + size_steps + 1;
        data(1:size_steps,1:3,i) = rawData(pos-(size_steps+1)+2:pos, 1:3);
    end
    
%     t = data(:,1,1);
%     [~, timestep] = min(abs(t-time));
    
    step_sorted_data = [dt', reshape(data(end,2:3,:),[2,length(dt)])'];
    [Q(j), dQ(j)] = textra(step_sorted_data, 1, 0.25);
end

if(Rodlike_opt)
    %for shear-rate xvals
    xvals = xvals*4*sigma^2;
    %for Viscosity data
    Q = Q/(4*sigma^2);
    dQ = dQ/(4*sigma^2);
    %for Psi1, Psi2 data
%     Q = Q/(16*sigma^4);
%     dQ = dQ/(16*sigma^4);
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
% axes1.YScale='log';
e1 = errorbar(xvals, Q, dQ, 'bd', ...
        'DisplayName','TEXTRA','LineWidth',2);
e1.MarkerFaceColor='b';
e1.MarkerSize=14;
dim = [0.18 0.18 0.7 0.7];
% str = ['t = ', num2str(time)];
% annotation('textbox',dim,'String',str,'FitBoxToText',...
%     'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
hold off

