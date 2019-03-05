%% Time Data from timestepdata.inp
% Always run this first
clear variables
dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
% Nrelax_times = 1100;
delay = dtData(1,4);
size_steps=Nrelax_times/delay;
inpdata = dlmread('inputparameters.inp', '');
sigma = inpdata(4);
alpha = inpdata(2);
%set Rodlike_opt to true for Rodlike units, false for Hookean
Rodlike_opt = true;
warning('off', 'MATLAB:nearlySingularMatrix')
xvals = dlmread('srvals.dat', '');
tmaxvals = dlmread('tvals.dat','');

%% Read a single file and plot time-series data at different timesteps

file = 'eta_sr0.005.dat';
rawData = dlmread(file, '');
Nrelax_times = tmaxvals(xvals==0.005);

data = nan(size_steps,3,length(dt));
pos = 0;
for i=1:length(dt)
    if dt(i)>=delay
        size_steps=Nrelax_times/dt(i);
    else
        size_steps=Nrelax_times/delay;
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
% errorbar(data(:,1,i), Q, dQ,...
%         'DisplayName','TEXTRA','LineWidth',2);
[h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',16,'Interpreter','latex','Box','off');
hold off

%% Read a list of files and plot data at a specific time
% xvals = [0.0447676953 0.2238384766 0.4476769531 0.8953539063 1.3430308594];

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
    clear data
    rawData = dlmread(files(j), '');
    
    Nrelax_times = tmaxvals(j);
%     data = nan(size_steps,3,length(dt));
    pos = 0;
    for i=1:length(dt)
        if dt(i)>delay
            size_steps=Nrelax_times/dt(i);
        else
            size_steps=Nrelax_times/delay;
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
axes1.YScale='log';
e1 = errorbar(xvals, Q, dQ, 'gd', ...
        'DisplayName','FF Simulations','LineWidth',2);
e1.MarkerFaceColor='g';
e1.MarkerSize=14;
% dim = [0.18 0.18 0.7 0.7];
% str = {['\delta Q_R = ', num2str(0.005)],...
%        ['H_R = ', num2str(400)],...
%        ['N_{traj} = 10^7'],...
%        ['h_R = 3/8']};
% annotation('textbox',dim,'String',str,'FitBoxToText',...
%     'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
hold off

writeFile = fopen('H10-7_dQ10-4_eta.dat', 'w');
fprintf(writeFile, '%e, %e, %e\n', [xvals; Q; dQ]);
fclose(writeFile);
