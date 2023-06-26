% ================================================================== %
% This script verifies the correctness of Dijkstra_with_early_stopping
% Specifically, we compute, for all datapoints in a medium sized data set,
% their k nearest neighbors in a power weighted shortest path distance. 
% We do this in two ways: using Dijkstra_with_early_stopping and using
% Dijkstra, and verify that the outputs are the same for both cases.
% Daniel Mckenzie
% June 2019
% ================================================================== %

clear, close all, clc
addpath(genpath('../Utilities'))
addpath(genpath('../Modified_Dijkstra'))

% ========================= Generate Data ======================= %
ysep= 1; len = 3; n0 = 100; noise_level = 0.14; ambient_dim = 100;
k = 3;
N = k*n0;
colours = ['r', 'b','g', 'm', 'c', 'y',];


[Points2D,Points] = Generate3Lines(ysep,len,n0,noise_level,ambient_dim);
perm = randperm(N);
[~,permInv] = sort(perm);
Data = Points(perm,:);
Points2D = Points2D(perm,:);
plot(Points2D(:,1),Points2D(:,2),'k*')


% =========== Find nearest neighbors using naive Dijkstra ============== %
p = 10; % the power weighting
K = 25; % number of nearest neighbors to consider

tic
Dists = squareform(pdist(Data));
WeightedDists = Dists.^p;
AdjMat = Dists > 0;
[costs,paths] = dijkstra(AdjMat,WeightedDists);
[NaivePathDist,NaiveKNN]  = mink(costs,K);
NaiveTime = toc;

% ========= Find nearest neighbors using Dijkstra_with_pruning ======= %
tic
[IDX,D] = knnsearch(Data,Data,'K',K,'NSMethod','kdtree');
D = D.^p;  % power weight the distances

FastPathDist = zeros(K,N);
FastKNN = zeros(K,N);

for i = 1:N
    [K_nearest_neighbors, Distances] = Dijkstra_with_early_stopping(IDX,D,i);
    FastPathDist(:,i) = Distances;
    FastKNN(:,i) = K_nearest_neighbors;
end
FastTime = toc;

% =========== check to see if both methods find same results =========== %

KNN_is_equal = isequal(NaiveKNN,FastKNN);
Difference_between_dists = norm(NaivePathDist - FastPathDist, 'fro')

if KNN_is_equal && Difference_between_dists <=1e-10
    disp(['Both Approaches yielded identical output, although Naive Dijkstra ran in ',...
        num2str(NaiveTime),' seconds while Dijkstra_with_pruning only required ',num2str(FastTime),' seconds.'])
else
    disp(['Dijkstra_with_pruning did not yield the same output as the naive Dijkstra approach'])
end

% ============ Illustrate some of these nearest neighbors ======== %
Inds = 1:N;
hold on
for i = 1:6
    RandPoint = datasample(Inds,1);
    NearestNeighbors = FastKNN(:,RandPoint);
    plot(Points2D(RandPoint,1),Points2D(RandPoint,2),strcat(colours(i),'o'))
    plot(Points2D(NearestNeighbors,1),Points2D(NearestNeighbors,2),strcat(colours(i),'*'))
    Inds = setdiff(Inds,NearestNeighbors);
end
    