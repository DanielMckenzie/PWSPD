function A = SimilarityMatrixZMP(Data,r)
%        A = SimilarityMatrixZMP(Data,r)
% Create a full, ZMP scaled adjacency matrix.
% Do not use for large data sets!
% Daniel Mckenzie
% 21 June 2019
%
% INPUT
% ====================================================
% Data .................... n-by-d data matrix. Data points stored as rows
% r .................... Local clustering parameters will be set using the
% r-th nearest neighbour
% 
% OUTPUT
% =========================================
% A ................... Weighted adjacency matrix.

[n,d] = size(Data);
Dists  = squareform(pdist(Data));

Scales = mink(Dists,r);
ScaleMat = Scales'*Scales;
A = exp(-(Dists.^2)./ScaleMat);

end