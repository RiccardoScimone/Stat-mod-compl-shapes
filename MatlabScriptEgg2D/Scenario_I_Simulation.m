%%%%% This simulate wider struts (Scenario I). Code is the same that
%%%%% In-control simulation with mu_s1 instead of mu_s0
clear all;
clc;
close all;
load('Egg_analytical_model.mat');
nsim = 5; % Elements to be simulated
mu_s1 = 0.055; % Shifted mean, TO BE SET 
dim = size(vx);
dim = dim(2);
npoints_max = 1000; %points per segment
npoints = ceil(lengths/max(lengths)*npoints_max);

%%%%% Simulation parameters as in the work %%%%%%
theta = 0.008;
sigma_bar = 0.004;
epsilon = 0.002;
sigma1 = sigma_bar - epsilon;
sigma2 = sigma_bar + epsilon;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')
sigmarvec = sigma1 + rand(nsim,1)*(sigma2 - sigma1);%%%% generate sigma_ij for each shape
parfor k = 1:nsim %%%% for each shape to be simulated
    widths = mu_s1 + normrnd(0,theta,dim,1); %%%% Generate the width of each strut within the shape
    PointCloud = [x y]; %%%% This will be the final shape. [x y] are the boundary poinits.
    for i = 1:dim %For each strut
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
%% Comment/uncomment for plotting and saving
   figure;
   scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
   axis equal;
  writematrix(PointCloud,strcat('PointClouds/OOC_WS',num2str(k),'.csv'));
end

