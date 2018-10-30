%% Setup
clear variables

dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
Ntraj = dtData(1,1);
delay = dtData(1,4);
size_steps=floor(Nrelax_times/delay+1);
fileID = fopen('Q_dist_output.dat', 'r');
data = ones(Ntraj, 3, length(dt), size_steps);

tic
for i=1:length(dt)
    dummy = fgets(fileID);
    if dt(i)>=delay
        size_steps=floor(Nrelax_times/dt(i));
    else
        size_steps=floor(Nrelax_times/delay+1);
    end
    
    for j=1:size_steps
        if (j~=1)
            dummy=fgets(fileID);
        end
        dummy = fgets(fileID);
        data(:,:,i,j) = fscanf(fileID, '%f %f %f', size(data(:,:,i,j)'))';
    end
end
toc

fclose(fileID);

%% Processing and plotting

for i=1:Ntraj
    x = data(i,1,1,1);
    y = data(i,2,1,1);
    z = data(i,3,1,1);
    [azi(i), elev(i), r(i)] = cart2sph(x,y,z);
end

% for some reason this isn't giving the exact correct results.
azi_range = -pi:pi/12:pi;
elev_range = -pi/2:pi/12:pi/2;
[azi_counts, azi_range] = histcounts(azi, azi_range, 'Normalization', 'pdf');

azi_mid = (azi_range(1:end-1)+azi_range(2:end))/2;
figure();
polarplot(azi_mid, azi_counts);

[counts, azi_r, elev_r] = histcounts2(azi, elev, azi_range, elev_range,...
    'Normalization', 'pdf');

azi_m = (azi_r(1:end-1)+azi_r(2:end))/2;
elev_m = (elev_r(1:end-1)+elev_r(2:end))/2;

[azi_g, elev_g] = meshgrid(azi_m, elev_m);

[X,Y,Z] = sph2cart(azi_g', elev_g', counts);

figure();
surf(X,Y,Z);
axis square
