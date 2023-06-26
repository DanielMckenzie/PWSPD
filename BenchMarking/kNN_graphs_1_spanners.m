% ====================================================================== %
% Experimentally verifying that the k-NN graph is a 1 spanner
% Part of code for:
%     " Geometry, Density and Path Distances on Graphs"
%               by Little, McKenzie and Murphy
% Daniel McKenzie
% March 02 2020
% ======================================================================= %

%clear, close all, clc
addpath(genpath('../Utilities'),genpath('../ClusteringAlgorithms'))
addpath(genpath('../Modified_Dijkstra'),genpath('../PathMetric'))
addpath('../MAT_Files')


% ======================= Generate the Data ============================= %
n = 100;
d = 5;
Data = rand(n,d); % uniformly distributed data on unit cube in R^d.
%Data = 0.5*randn(n,d);
% Step-up distribution
%pd = makedist('PieceWiseLinear',[0 0.5 1], [0 0.2 1]);
%x1 = random(pd,n,1);
%Data = [x1, rand(n,d-1)];
p = 1;


% =============== Create the k nearest neighbor graphs ================= %
eps = (1/(4^(1/p)) - 1/4)^(d/2)
%k = ceil(1-3*log(n)/log(1-eps))
k = 100
A_k = CreateKNN_Max_from_Data2(Data,k).^p;
AdjMat = A_k > 0;
[k_costs,k_paths] = dijkstra(AdjMat,A_k);
[k_costs_vec, sort_order] = sort(reshape(k_costs,n^2,1),'descend');
figure
plot(1:n^2,k_costs_vec,'b')
hold on

% =================== Create the full graphs =========================== %
Dists = squareform(pdist(Data));
WeightedDists = Dists.^p;
AdjMat = Dists > 0;
[full_costs,full_paths] = dijkstra(AdjMat,WeightedDists);
full_costs_vec = reshape(full_costs,n^2,1);
full_costs_vec = full_costs_vec(sort_order);
%figure
plot(1:n^2,full_costs_vec,'r--')
legend('Path distance in k-NN graph','Path distance in complete graph')
set(gca,'FontSize',15)

rel_error = (k_costs_vec - full_costs_vec)./full_costs_vec;
figure
plot(1:n^2,rel_error)
