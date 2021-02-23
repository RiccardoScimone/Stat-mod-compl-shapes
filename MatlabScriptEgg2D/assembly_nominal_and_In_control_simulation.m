%%% This script gemerates nsim "in control" eggs.
%%%% nsim is always 95 in the work. The name of the variables
%%%% are coherent with the paper.
%%%% The shapes are symulated from the analytical model in 
%%%% Egg_analytical_model.mat
clear all;
clc;
close all;
load('Egg_analytical_model.mat'); % Informations about the analytical model 
dim = size(vx); % Number of segments. vx and vy are the coordinate of the endpoints of each vertex
dim = dim(2); %number of segments, vx and vy are the vertices coordinates
npoints_max = 1000; %points per segment
npoints = ceil(lengths/max(lengths)*npoints_max);
nsim = 95; %in control elements

%%%%% Simulation parameters as in the work %%%%%%%
mu_s0 = 0.04;
theta = 0.008;
sigma_bar = 0.004;
epsilon = 0.002;
sigma1 = sigma_bar - epsilon;
sigma2 = sigma_bar + epsilon; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle');
sigmarvec = sigma1 + rand(nsim,1)*(sigma2 - sigma1); %%%% generate sigma_ij for each shape
%% Assembling In Control Elements
parfor k = 1:nsim %drop parallelization if you want plots
    widths = mu_s0 + normrnd(0,theta,dim,1); %%%% Generate the width of each strut within the shape
    PointCloud = [x y]; %%%% This will be the final shape. [x y] are the boundary poinits.
    for i = 1:dim %create struts
        phi = atan((vx(2,i)-vx(1,i))/(vy(2,i)-vy(1,i))); 
        rot = [cos(phi) -sin(phi); sin(phi) cos(phi)]; 
        invrot = rot'; 
        vert = [vx(2,i)-vx(1,i) vy(2,i)-vy(1,i)]';
        vert = rot*vert;
        vert(1) = 0;
        vdisc = [(zeros(npoints(i),1)) (linspace(0,vert(2),npoints(i)))'];
        for j = 1:npoints(i)
            vdisc(j,1) = ((-1)^j)*widths(i)/2+normrnd(0,sigmarvec(k));
        end
        vdisc = (invrot*vdisc')';
        for j = 1:npoints(i)
            vdisc(j,:) = vdisc(j,:) + [vx(1,i) vy(1,i)];
        end
        PointCloud = [PointCloud; vdisc];
    end
    for j = 1:ninv
        PointCloud(j,:) = PointCloud(j,:) - abs(normrnd(0,sigma_bar))*neggfun(T(j))';
    end
%% Comment/Uncomment for saving or plotting
%    figure;
%    scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
%    axis equal;
writematrix(PointCloud,strcat('PointClouds/In_control_',num2str(k),'.csv'));
end
%% Assembling nominal model
PointCloud = [x y]; 
for i = 1:dim 
    phi = atan((vx(2,i)-vx(1,i))/(vy(2,i)-vy(1,i))); 
    rot = [cos(phi) -sin(phi); sin(phi) cos(phi)]; 
    invrot = rot';
    vert = [vx(2,i)-vx(1,i) vy(2,i)-vy(1,i)]';
    vert = rot*vert;
    vert(1) = 0;
    vdisc = [(zeros(npoints(i),1)) (linspace(0,vert(2),npoints(i)))'];
    for j = 1:npoints(i)
        vdisc(j,1) = ((-1)^j)*mu_s0/2+normrnd(0,sigma_bar);
    end
    vdisc = (invrot*vdisc')';
    for j = 1:npoints(i)
        vdisc(j,:) = vdisc(j,:) + [vx(1,i) vy(1,i)];
    end
    PointCloud = [PointCloud; vdisc];
end
for j = 1:ninv
    PointCloud(j,:) = PointCloud(j,:) - abs(normrnd(0,sigma_bar))*neggfun(T(j))';
end
%% Comment/uncomment for saving or plotting
%     figure;
%     scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
%     xlim([min(x) max(x)]);
%     ylim([min(y) max(y)]);
%     axis equal;
%     title('Point Cloud Example')
     writematrix(PointCloud,strcat('PointClouds/Nominal_model.csv'));
