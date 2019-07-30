%% Time Data from timestepdata.inp
% Always run this first
clear variables
dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
Ntraj = dtData(1,1);
% Nrelax_times = 1100;
delay = dtData(1,4);
size_steps=round(Nrelax_times/delay);
inpdata = dlmread('inputparameters.inp', '');
sigma = inpdata(4);
alpha = inpdata(2);
%set Rodlike_opt to true for Rodlike units, false for Hookean
Rodlike_opt = true;
warning('off', 'MATLAB:nearlySingularMatrix')
% xvals = dlmread('srvals.dat', '');
xvals=[0 0.02 0.05 0.1 0.15 0.25 0.375];
% tmaxvals = dlmread('tvals.dat','');
tmaxvals=ones(size(xvals))*80;

%% Read a list of files and plot data at a specific time
% xvals = [0.0447676953 0.2238384766 0.4476769531 0.8953539063 1.3430308594];

% If you've used a different naming convention, you can just set the 
% files variable as an array of strings directly
prefix = "eta_hstar";
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

yvals = [];
dyvals = [];
clear Q dQ

for j = 1:length(files)
    file = files(j);
    rawData = dlmread(file, '');
    Nrelax_times = tmaxvals(j);
    size_steps=round(Nrelax_times/delay);

    data = nan(size_steps,3,length(dt));
    pos = 0;
    for i=1:length(dt)
        if dt(i)>=delay
            size_steps=Nrelax_times/dt(i);
        else
            size_steps=round(Nrelax_times/delay);
        end
        pos = pos + size_steps + 1;
        data(1:size_steps,1:3,i) = rawData(pos-(size_steps+1)+2:pos, 1:3);
    end

    if(Rodlike_opt)
        data(:,1,:) = data(:,1,:)/(4*sigma^2);
        %for Viscosity data
        data(:,2:3,:) = data(:,2:3,:)/(4*sigma^2);
%         Q = Q/(4*sigma^2);
%         dQ = dQ/(4*sigma^2);
        %for Psi1, Psi2 data
    %     data(:,2:3,:) = data(:,2:3,:)/(16*sigma^4);
    %     Q = Q/(16*sigma^4);
    %     dQ = dQ/(16*sigma^4);

    end

    fig = figure();
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
    
    [t,val] = ginput(1);
    
    %actually won't work if size_steps changes, aka if dt(i)>delay
    no_cells = size_steps - nnz(data(:,1,1)>t);
    
    for i = 1:length(dt)
        avg(j,i) = mean(data(no_cells:end,2,i));
        vars = (data(no_cells:end,3,i)*sqrt(Ntraj)).^2;
        devs = (data(no_cells:end,3,i) - avg(j,i)).^2;
        total_var = sum(vars + devs)/(size_steps-no_cells);
        err(j,i) = sqrt(total_var/(Ntraj*(size_steps-no_cells)));
%         err_check(j,i) = std(data(no_cells:end,3,i))/sqrt(size_steps-no_cells);
    end
    
    [Q(j), dQ(j)] = textra([dt',avg(j,:)',err(j,:)'], 1, 0.1);
    
    close(fig)
    
    clear dt_data
    
end

% figure();
% hold on
% axes1 = gca;
% fsize=20;
% pbaspect([1. 1. 1]);
% hold(axes1,'on');
% box(axes1,'on');
% set(axes1,'FontSize',16,'LineWidth',2,'TickLength',[0.015 0.025]);
% % Comment out the below 2 lines
% % as needed.
% axes1.XScale='log';
% axes1.YScale='log';
% e1 = errorbar(xvals, Q, dQ, 'bd', ...
%         'DisplayName','$H^*_R = 10, s=0.01$','LineWidth',2);
% e1.MarkerFaceColor='b';
% e1.MarkerSize=14;
% % dim = [0.18 0.18 0.7 0.7];
% % str = {['\delta Q_R = ', num2str(0.005)],...
% %        ['H_R = ', num2str(400)],...
% %        ['N_{traj} = 10^7'],...
% %        ['h_R = 3/8']};
% % annotation('textbox',dim,'String',str,'FitBoxToText',...
% %     'on','Interpreter','latex','FontSize',14,'EdgeColor','None');
% hold off

% writeFile = fopen('H10_s0.01_chiTau.dat', 'w');
% fprintf(writeFile, '%e, %e, %e\n', [xvals; Q; dQ]);
% fclose(writeFile);

figure();
hold on
for i=1:length(dt)
    name = ['dt=',num2str(dt(i))];
    e1 = errorbar(xvals, avg(:,i), err(:,i), 'd', 'DisplayName', name);
end
e1 = errorbar(xvals, Q, dQ, 'd', 'DisplayName', 'TEXTRA');
[h,icons,plots,legend_text]=legend({},'Location','northwest','FontSize',16,'Interpreter','latex','Box','off');
hold off

