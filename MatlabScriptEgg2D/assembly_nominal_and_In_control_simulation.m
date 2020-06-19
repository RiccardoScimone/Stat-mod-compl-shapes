%%% This script gemerates nsim "in control" eggs.
%%%% nsim is always 95 in the work. The name of the variables
%%%% are coherent with the paper.
%%%% The shapes are symulated from the analytical model in 
%%%% Egg_analytical_model.mat
clc;
clear all;
close all;
load('Egg_analytical_model.mat'); % Informations about the analytical model 
dim = size(vx); 
dim = dim(2); % Number of segments. vx and vy are the coordinate of the endpoints of each vertex
npoints = 400; %number of points which discretizes each segment. 
nsim = 95; % number of simulated in control shapes
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

parfor k = 1:nsim %%%% for each shape to be simulated
widths = mu_s0 + normrnd(0,theta,dim,1); %%%% Generate the width of each strut within the shape
PointCloud = [x y]; %%%% This will be the final shape. [x y] are the boundary poinits.
for i = 1:dim %For each strut
    phi = atan((vx(2,i)-vx(1,i))/(vy(2,i)-vy(1,i))); %%angle of the strut wrt y axis
    rot = [cos(phi) -sin(phi); sin(phi) cos(phi)]; %%corresponding rotation matrix
    invrot = rot'; %%inverse rotation
    %center and rotate the strut on y axis
    vert = [vx(2,i)-vx(1,i) vy(2,i)-vy(1,i)]'; 
    vert = rot*vert;
    vert(1) = 0;
    %%%%%
    %discretize the vertical segment
    vdisc = [(zeros(npoints,1)) (linspace(0,vert(2),npoints))'];
    %%%%
    %double the segment creating the strut, apply perturbation
     for j = 1:npoints
         vdisc(j,1) = ((-1)^j)*widths(i)/2+normrnd(0,sigmarvec(k));
     end
     %%%%%
     %back rotate and translate
     vdisc = (invrot*vdisc')';
     for j = 1:npoints
         vdisc(j,:) = vdisc(j,:) + [vx(1,i) vy(1,i)];
     end
     %%%%%%
     %%% put the strut in the point cloud
     PointCloud = [PointCloud; vdisc];
     %%%%%%%
end
%Perturbation of the boundary
for j = 1:ninv
    PointCloud(j,:) = PointCloud(j,:) - abs(normrnd(0,sigma_bar))*neggfun(T(j))';
end
%%%% uncomment if you want the plots of each shape. It must stay commented
%%%% when running from R.
  %  figure;
  %  scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
  %  axis equal;
%%%%% Save the point cloud
    writematrix(PointCloud,strcat('PointClouds/In_control',num2str(k),'.csv'));
end



%%%%% assembly nominal model with the same procedure 
 PointCloud = [x y]; 
 for i = 1:dim 
    phi = atan((vx(2,i)-vx(1,i))/(vy(2,i)-vy(1,i))); 
    rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    invrot = rot'; 
    vert = [vx(2,i)-vx(1,i) vy(2,i)-vy(1,i)]';
    vert = rot*vert;
    vert(1) = 0;
    vdisc = [(zeros(npoints,1)) (linspace(0,vert(2),npoints))'];
     for j = 1:npoints
         vdisc(j,1) = ((-1)^j)*mu_s0/2+normrnd(0,sigma_bar);
     end
     vdisc = (invrot*vdisc')';
     for j = 1:npoints
         vdisc(j,:) = vdisc(j,:) + [vx(1,i) vy(1,i)];
     end
     PointCloud = [PointCloud; vdisc];
 end
for j = 1:ninv
    PointCloud(j,:) = PointCloud(j,:) - abs(normrnd(0,sigma_bar))*neggfun(T(j))';
end
     %figure;
     %scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
     %xlim([min(x) max(x)]);
     %ylim([min(y) max(y)]);
     %axis equal;
     %title('Point Cloud Example')
    writematrix(PointCloud,strcat('PointClouds/Nominal_model.csv'));
