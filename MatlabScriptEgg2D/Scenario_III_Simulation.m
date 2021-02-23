%% This script generates Out of Control elements, with irregularities
% Parameters meaning is the same as in control simulations
clear all;
clc;
close all;
load('Egg_analytical_model.mat');
dim = size(vx);
dim = dim(2);
npoints_max = 1000; %points per segment
npoints = ceil(lengths/max(lengths)*npoints_max);
nsim = 5;
mu_s0 = 0.04;
theta = 0.008;
sigma_bar = 0.004;
epsilon = 0.002;
sigma1 = sigma_bar - epsilon;
sigma2 = sigma_bar + epsilon; 
sigmarvec = sigma1 + rand(nsim,1)*(sigma2 - sigma1);
irrsigma = 0.001; %%%% a small sigma to make the irregularities bidimensional
%% This tunes irregularities
nu = 10; %%TO BE SET
minfrac = 1/20;
maxfrac = 1/16;
%fil_info = randi([nfil/2,nfil],1,length(lengths));
%placefilm = rand(nsim,sum(fil_info));
%% Now the work is as before
parfor k = 1:nsim
    fil_info = randi([0,nu],1,length(lengths));
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
        realnfil = fil_info(i);% Quanti filamenti?
        placefil = rand(1,fil_info(i));
        seglength = sqrt((vx(2,i)-vx(1,i))^2+(vy(2,i)-vy(1,i))^2);
        for r = 1:realnfil %Ogni filamento deve essere piazzato
            fillength = (minfrac + rand*(maxfrac - minfrac))*seglength; %lunghezza del filamento
            npointsfil = floor(npoints(i)*fillength/seglength/2);
            yfil = vy(1,i)+placefil(r)*(vy(2,i)-vy(1,i));
            xfil = vx(1,i)+placefil(r)*(vx(2,i)-vx(1,i));
            v = [xfil - vx(1,i);yfil-vy(1,i)];
            vrot = rot*v;
            temp = rand;
            fildisc = [((temp > 0.5)*2 - 1)*widths(i)/2+((temp > 0.5)*2 - 1)*linspace(0,fillength,npointsfil)' vrot(2)*ones(npointsfil,1)+normrnd(0,irrsigma,npointsfil,1)];
            fildisc = (invrot*fildisc')';
            for j = 1:npointsfil
                fildisc(j,:) = fildisc(j,:) + [vx(1,i) vy(1,i)];
            end
            PointCloud = [PointCloud; fildisc]; 
        end   
    end
    %Perturbazione dell'involucro, VERSO L'ESTERNO
    for j = 1:ninv
        PointCloud(j,:) = PointCloud(j,:) - abs(normrnd(0,sigma_bar))*neggfun(T(j))';
    end
    %% Comment/uncomment for saving or plotting
%     figure;
%     scatter(PointCloud(:,1), PointCloud(:,2),1,'filled');
%     xlim([min(x) max(x)]);
%     ylim([min(y) max(y)]);
%     axis equal;
    writematrix(PointCloud,strcat('PointClouds/OOC_IRR',num2str(k),'.csv'));
end

