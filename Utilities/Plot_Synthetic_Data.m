% =================================================================== %
% Simple script for plotting synthetic data sets.
% Daniel Mckenzie
% 9th September 2019
% =================================================================== %

clear, close all, clc

% =========================== Parameters ========================= %
colours = ['r', 'b','g','k', 'm', 'c', 'y',];
n0 = 500; noise_level = 0.14; k=3; ambient_dim = 50;
num_trials = 10;
sym = 'max';

% == Moons
% r1 = 1;
% r2 = 1;
% r3 = 1.5;

% == Lines
%ysep=1; len = 5;

% == Circles
r1 = 1;
r2 = 2.25;
r3 = 3.5;
n1 =222;
n2 = 500;
n3 = 778;
N = n1+n2+n3;

%N = k*n0;

% == Lines
%[Points2D,Points] = Generate3Lines(ysep,len,n0,noise_level,ambient_dim);
% == Moons
%[Points2D,Points] = Generate3Moons(r1,r2,r3,n0,noise_level,ambient_dim);
% == Circles
[Points2D,Points] = Generate3Circles(r1,r2,r3,n1,n2,n3,noise_level,ambient_dim);

plot(Points2D(:,1),Points2D(:,2), 'b*')
set(gca,'FontSize',18)

    
    