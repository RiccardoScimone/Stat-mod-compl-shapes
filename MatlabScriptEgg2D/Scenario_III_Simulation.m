%%% This script generates eggs with irregularities as in scenario 3
%%% All parameters are named as in the work
clc;
clear all;
close all;
nu = 20; %Max number of irregularities, TO BE SET
load('Egg_analytical_model.mat');
dim = size(vx);
dim = dim(2); 
npoints = 400; 
nsim = 5;
mu_s0 = 0.04;
theta = 0.008;
sigma_bar = 0.004;
epsilon = 0.002;
sigma1 = sigma_bar - epsilon;
sigma2 = sigma_bar + epsilon; 
sigmarvec = sigma1 + rand(nsim,1)*(sigma2 - sigma1);
irrsigma = 0.001; %%%% a small sigma to make the irregularities bidimensional
rng('shuffle');
for i = 1:dim
   lengths(i) = sqrt((vx(2,i)-vx(1,i))^2+(vy(2,i)-vy(1,i))^2); %%% Length of each strut
end
minfrac = 1/10; %minimum fraction of the length of the strut where the irregularity is attached
maxfrac = 1/8; %maximum fraction of the length of the strut where the irregularity is attached
parfor k = 1:nsim 
irr_info = randi([nu/2,nu],1,dim); %% irregularities for each strut 
widths = mu_s0 + normrnd(0,theta,dim,1);
PointCloud = [x y]; 
for i = 1:dim %per ogni segmento
    placeirrs = rand(irr_info(i));% Tell where each irregularity has to be placed
    phi = atan((vx(2,i)-vx(1,i))/(vy(2,i)-vy(1,i))); 
    rot = [cos(phi) -sin(phi); sin(phi) cos(phi)]; 
    invrot = rot'; 
    vert = [vx(2,i)-vx(1,i) vy(2,i)-vy(1,i)]';
    vert = rot*vert;
    vert(1) = 0;
    vdisc = [(zeros(npoints,1)) (linspace(0,vert(2),npoints))'];
     for j = 1:npoints
         vdisc(j,1) = ((-1)^j)*widths(i)/2+normrnd(0,sigmarvec(k));
     end
     vdisc = (invrot*vdisc')';
     for j = 1:npoints
         vdisc(j,:) = vdisc(j,:) + [vx(1,i) vy(1,i)];
     end
     PointCloud = [PointCloud; vdisc];
     %%% now the irregularities
     nirrs = irr_info(i);
     for r = 1:nirrs %%% For each irregularity
        irrlength = (minfrac + rand*(maxfrac - minfrac))*lengths(i); %length of the irregularity
        npointsirr = floor(npoints*irrlength/lengths(i)); % number of points discretizing the irregularity
        placeirr = placeirrs(r); %%%where to place the irregularity
        %%%% now consruct the irregularity
        yirr = vy(1,i)+placeirr*(vy(2,i)-vy(1,i));
        xirr = vx(1,i)+placeirr*(vx(2,i)-vx(1,i));
        v = [xirr - vx(1,i);yirr-vy(1,i)];
        vrot = rot*v;
        temp = rand;
        irrdisc = [((temp > 0.5)*2 - 1)*widths(i)/2+((temp > 0.5)*2 - 1)*linspace(0,irrlength,npointsirr)' vrot(2)*ones(npointsirr,1)+normrnd(0,irrsigma,npointsirr,1)];
        irrdisc = (invrot*irrdisc')'; 
        for j = 1:npointsirr
         irrdisc(j,:) = irrdisc(j,:) + [vx(1,i) vy(1,i)];
        end
        %%%%%%
        PointCloud = [PointCloud; irrdisc]; %Append the irregularity to the point cloud        
    end
end
%%%% boundary of the shape
for j = 1:ninv
    PointCloud(j,:) = PointCloud(j,:) - abs(normrnd(0,sigma_bar))*neggfun(T(j))';
end
%figure;
%scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
%xlim([min(x) max(x)]);
%ylim([min(y) max(y)]);
%axis equal;
writematrix(PointCloud,strcat('PointClouds/OOC_IRR',num2str(k),'.csv'));
end

