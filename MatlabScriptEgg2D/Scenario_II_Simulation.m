%Simulate nsim out of controls with nu missing struts (Scenario II)
% Code is exactly the same that the In-control simulation,
% We just eliminate the struts first
clc;
clear all;
close all;
eta = 4; %%% number of missing trabecula, TO BE SET
load('Egg_analytical_model.mat');
dim = size(vx);
dim = dim(2);
npoints_max = 1000; %points per segment
npointsconst = ceil(lengths/max(lengths)*npoints_max);
nsim = 5;

mu_s0 = 0.04;
theta = 0.008;
sigma_bar = 0.004;
epsilon = 0.002;
sigma1 = sigma_bar - epsilon;
sigma2 = sigma_bar + epsilon;
rng('shuffle');
sigmarvec = sigma1 + rand(nsim,1)*(sigma2 - sigma1);
vxcost = vx;
vycost = vy;
dimcost = dim;

parfor k = 1:nsim
    vx = vxcost; %drop struts
    vy = vycost;
    npoints = npointsconst;
    dim = dimcost;
    todelete = randi([1 dim],1,eta);
    vx(:,todelete) = [];
    vy(:,todelete) = [];
    npoints(todelete) = [];
    dim = size(vx);
    dim = dim(2); %struts dropped 
    widths = mu_s0 + normrnd(0,theta,dim,1);
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
%%  Comment/uncomment for plotting and saving
%   figure;
%   scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
%   axis equal;
%   title('4 missing')
  writematrix(PointCloud,strcat('PointClouds/OOC_MS',num2str(k),'.csv'));
end
