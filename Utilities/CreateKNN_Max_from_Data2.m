function A = CreateKNN_Max_from_Data2(Data,k)
%        A = CreateKNN_Max_from_Data2(Data,k)
% Create the adjacency matrix of a (weighted) K-Nearest Neighbours graph
% from a DATA MATRIX X. Gaussian kernel is used.
% No local scaling.
% Symmetrize by taking a max.
% 
% INPUT
% =========================================
% Data ................. n-by-d data matrix. Data points stored as rows
% K .................... Number of nearest neighbours to keep.
%
%OUTPUT
% =========================================
% A ................... Weighted adjacency matrix.
%
% Daniel Mckenzie, using some code from Ke Yin's graph_affinity_matrix.m
% 18th December 2018

[n,d] = size(Data);
[IDX,D] = knnsearch(Data,Data,'K',k,'NSMethod','kdtree');
Dists = reshape(D',[],1);
I = reshape(ones(k, 1) * (1:n), [], 1);
J = reshape(IDX',[],1);

Atemp = sparse(I,J,Dists,n,n);
A = max(Atemp',Atemp);

end

    
    