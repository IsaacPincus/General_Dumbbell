dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
delay = dtData(1,4);
size_steps=Nrelax_times/delay+1;

rawData = dlmread('eta.dat', '');

data = [];
for i=1:length(dt)
    data(:,:,i) = rawData((size_steps+1)*(i-1)+2:(size_steps+1)*i, 1:3);
end

step_sorted_data = [dt', reshape(data(11,2:3,:),[2,length(dt)])'];
textra_f_plot_data(step_sorted_data, 1, 0.1);

for i=1:size_steps
    step_sorted_data = [dt', reshape(data(i,2:3,:),[2,length(dt)])'];
    [Q(i), dQ(i)] = textra_f_data(step_sorted_data, 1, 0);
end

figure();
hold on
for i=1:length(dt)
    errorbar(data(:,1,i), data(:,2,i), data(:,3,i))
end
errorbar(data(:,1,i), Q, dQ)
legend(['dt=',num2str(dt(1))],...
       ['dt=',num2str(dt(2))],...
       ['dt=',num2str(dt(3))],...
       'TEXTRA');
hold off

