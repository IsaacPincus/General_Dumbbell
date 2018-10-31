%% Setup
clear variables

dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
Ntraj = dtData(1,1);
delay = dtData(1,4);
size_steps=floor(Nrelax_times/delay+1)+1;
fileID = fopen('Q_dist_output.dat', 'r');
data = ones(Ntraj, 3, length(dt), size_steps);

tic
for i=1:length(dt)
    dummy = fgets(fileID);
    if dt(i)>=delay
        size_steps=floor(Nrelax_times/dt(i))+1;
    else
        size_steps=floor(Nrelax_times/delay+1)+1;
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

k = [0 1 0; 0 0 0; 0 0 0];
minX = -0.2;
maxX = 0.2;
minY = -0.2;
maxY = 0.2;
minZ = -0.2;
maxZ = 0.2;

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

tic
[X,Y,Z] = getSurf(data, 1, 1, Ntraj);
toc

% 3D figure with shear flow lines
aniFig = figure('units','normalized','outerposition',[0 0 1 1]);
axis equal
hold on
xlabel('x')
ylabel('y')
zlabel('z')
view(30,30)
hAni = surf(X,Y,Z);
% surf(X,Y,Z);
quiver3(Xs,Ys,Zs,...
    squeeze(fl(1,:,:,:)),squeeze(fl(2,:,:,:)),squeeze(fl(3,:,:,:)),...
    'linewidth', 2, 'MaxHeadSize', 0.5);

Frames(1) = getframe(aniFig);

pause(3)
for time = 2:size_steps
    [X,Y,Z] = getSurf(data, 1, time, Ntraj);
    hAni.XData = X;
    hAni.YData = Y;
    hAni.ZData = Z;
    drawnow limitrate
    Frames(time) = getframe(aniFig);
    
end
drawnow
hold off

movie2gif(Frames,'dist_evo.gif','DelayTime', 0.1)

%% Functions

function [X, Y, Z] = getSurf(timeData, timestep, time, Ntraj)
%     figure();
%     hold on
    for i=1:Ntraj
        x = timeData(i,1,timestep,time);
        y = timeData(i,2,timestep,time);
        z = timeData(i,3,timestep,time);
        [azi(i), elev(i), r(i)] = cart2sph(x,y,z);
%         scatter3(x,y,z)
    end
%     hold off
    
    
    
    azi_range = -pi:pi/12:pi;
    elev_range = -pi/2:pi/12:pi/2;

    [counts, azi_r, elev_r] = histcounts2(azi, elev, azi_range, elev_range);%,...
%     'Normalization', 'pdf')

    azi_m = (azi_r(1:end-1)+azi_r(2:end))/2;
    elev_m = (elev_r(1:end-1)+elev_r(2:end))/2;

    [azi_g, elev_g] = meshgrid(azi_m, elev_m);
    
    %Scale counts to account for equiangle grid bins
    counts = counts./abs(cos(elev_m));
    counts = counts/(sum(sum(counts))*(azi_r(1)-azi_r(2))*(elev_r(1)-elev_r(2)));
%     counts = counts/sum(sum(counts));

    [X,Y,Z] = sph2cart(azi_g', elev_g', counts);
    X = [X;X(1,:)];
    Y = [Y;Y(1,:)];
    Z = [Z;Z(1,:)];
end

