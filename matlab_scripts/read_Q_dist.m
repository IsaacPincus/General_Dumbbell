%% Setup
clear variables

dtData = dlmread('timestepdata.inp', '');
dt = dtData(2:end, 1)';
Nrelax_times = dtData(1, 3);
Ntraj = dtData(1,1);
delay = dtData(1,5);
size_steps=floor(Nrelax_times/delay+1);
fileID = fopen('Q_dist_output.dat', 'r');

inpdata = dlmread('inputparameters.inp', '');
sigma = inpdata(4);
alpha = inpdata(2);

tic
data = nan(Ntraj, 3, length(dt), size_steps);

file = 'Q_dist_output_sr0.dat';
fid = fopen(file, 'rb');

for i=1:length(dt)
    % read the first timestep width
    fprintf('Reading dt = %d \n', dt(i));
    hr1 = fread(fid, 1, 'int32');
    dt_check(i) = fread(fid, 1, 'float64');
    hr1 = fread(fid, 1, 'int32');
    
    % Find the length between outputs
    if dt(i)>=delay
        size_steps=Nrelax_times/dt(i)+1;
    else
        size_steps=Nrelax_times/delay+2;
    end
    
    % Run loop over timesteps
    for j=1:size_steps
        
        
        % read the first timestep
        hr1 = fread(fid, 1, 'int32');
        t(i,j) = fread(fid, 1, 'float64');
        hr1 = fread(fid, 1, 'int32');
        
        if (mod(j,10)==1)
            fprintf('reading t = %d \n', t(i,j));
        end

        % read the first set of data
        hrd = fread(fid, 1, 'int32');
        len = hrd/(3*8);
        data(:,:,i,j) = fread(fid, [3, len], 'float64')';
        hrd = fread(fid, 1, 'int32');
    end
end

fclose(fid);
toc

%% Histogram plot

width = 2*sqrt(alpha)/(100-1);
Q = max(0,(sigma-sqrt(alpha))):width:(sigma+sqrt(alpha));
Q(1) = Q(1) + 0.00000001;
Q(end) = Q(end) - 0.00000001;

phiQ = -alpha/2*log(1-(Q-sigma).^2/alpha);
Jeq = trapz(Q, exp(-phiQ).*Q.^2);
% Jeq = (1/(alpha+3)+sigma.^2/alpha)*beta(1/2,(alpha+2)/2)*alpha^1.5;
psiQ = Q.^2.*exp(-phiQ)./Jeq;

i = 4;

aniFigDist = figure();
hold on
plot(Q, psiQ, 'k-', 'Linewidth', 2, 'DisplayName', 'Analytical PDF at EQ');
j = 5;
Ql = sqrt(data(:,1,i,j).^2 + data(:,2,i,j).^2 + data(:,3,i,j).^2);
h1 = histogram(Ql, Q, 'Normalization', 'pdf');

FramesDist(1) = getframe(aniFigDist);

for j = 2:size_steps-1
    Ql = sqrt(data(:,1,i,j).^2 + data(:,2,i,j).^2 + data(:,3,i,j).^2);
    h1.Data = Ql;
    drawnow limitrate
    FramesDist(j) = getframe(aniFigDist);
end

hold off

movie2gif(FramesDist,'Q_dist_evo.gif','DelayTime', 0.5)

%% testing F(X) results

neg_fx = 0;

for i = 1:1
    for j=1:size_steps-1
        Ql = sqrt(data(:,1,i,j).^2 + data(:,2,i,j).^2 + data(:,3,i,j).^2);
        Fx = Ql + dt(i)/2*((Ql-sigma)./(1-(Ql-sigma).^2)/alpha);
        neg_fx = neg_fx + sum(Fx(Fx<0));
    end
end

%% Processing and plotting

timestep = 4;

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
[X,Y,Z] = getSurf(data, timestep, 1, Ntraj);
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
for time = 2:size_steps-1
    [X,Y,Z] = getSurf(data, timestep, time, Ntraj);
    hAni.XData = X;
    hAni.YData = Y;
    hAni.ZData = Z;
    drawnow limitrate
    Frames(time) = getframe(aniFig);
    pause(1)
    
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

