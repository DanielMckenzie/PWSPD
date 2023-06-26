% ====================================================================== %
% Benchmarking spectral clustering with Euclidean k-nn adjacency matrix vs.
% path distance k-nn adjacency matrix.
% Adjustable to all synthetic data sets
% Daniel Mckenzie
% 24th January 2019
% ======================================================================= %

clear, close all, clc
addpath(genpath('../Utilities'),genpath('../ClusteringAlgorithms'))
addpath(genpath('../Modified_Dijkstra'),genpath('../PathMetric'))
addpath('../MAT_Files')

% =========================== Parameters ========================= %
colours = ['r', 'b','g','k', 'm', 'c', 'y',];
n0 = 500; noise_level = 0.14; k=3; ambient_dim = 50;
num_trials = 10;
sym = 'max';

% == Moons
r1 = 1;
r2 = 1;
r3 = 1.5;

% == Lines
%ysep=1; len = 5;

% == Circles
% r1 = 1;
% r2 = 2.25;
% r3 = 3.5;
% n1 =222;
% n2 = 500;
% n3 = 778;
% N = n1+n2+n3;

N = k*n0;


% =================== Define all vectors of interest ============= %
Acc_Spectral_Path2_overall = 0;
Acc_Spectral_Path10_overall = 0;
Acc_Spectral_LLPD_overall = 0;
Acc_Spectral_Euclid_overall = 0;
Acc_Full_Spectral_Euclid_overall = 0;

Time_Spectral_Path2_overall = 0;
Time_Spectral_Path10_overall = 0;
Time_Spectral_LLPD_overall = 0;
Time_Spectral_Euclid_overall = 0;
Time_Full_Spectral_Euclid_overall = 0;


for i = 1:num_trials
    % ==================== Generate and Permute the data ============== %
    % == Lines
    %[Points2D,Points] = Generate3Lines(ysep,len,n0,noise_level,ambient_dim);
    % == Moons
    [Points2D,Points] = Generate3Moons(r1,r2,r3,n0,noise_level,ambient_dim);
    % == Circles
    %[Points2D,Points] = Generate3Circles(r1,r2,r3,n1,n2,n3,noise_level,ambient_dim);
    
    perm = randperm(N);
    [~,permInv] = sort(perm);
    Data = Points(perm,:);
    
    % =================== Find the ground truth clusters ============== %
    TrueClusters = cell(k,1);
    
    for a = 1:k
        TrueClusters{a} = permInv((a-1)*n0+1:a*n0);
    end
    
    % == Circles
    %TrueClusters{1} = permInv(1:n1);
  	%TrueClusters{2} = permInv(n1+1:n1+n2);
    %TrueClusters{3} = permInv(n1+n2+1:N);
    
     % ============== Full matrix Euclidean clustering =========== %
    tic
    A0 = SimilarityMatrixZMP(Data,10);
    [C0, ~, ~] = SpectralClustering(A0, k, 3);
    Time_Full_Spectral_Euclid_overall = Time_Full_Spectral_Euclid_overall + toc;
    ErrFullEuclid = 0;
    for a=1:k
        TempClust = find(C0(:,a));
        TempError = N;
        for b=1:k
            TempTempError = length(setdiff(TempClust,TrueClusters{b}));
            TempError = min(TempError,TempTempError);
        end
        ErrFullEuclid = ErrFullEuclid + TempError;
    end
    Acc_Full_Euclid_Cluster = 1 - ErrFullEuclid/N 
    
     % == Update
    Acc_Full_Spectral_Euclid_overall = Acc_Full_Spectral_Euclid_overall + Acc_Full_Euclid_Cluster;
   
    % ============== Euclidean knn clustering =========== %
    tic
    A1 = CreateKNN_Max_from_Data(Data,15,10);
    [C1, ~, ~] = SpectralClustering(A1, k, 3);
    Time_Spectral_Euclid_overall = Time_Spectral_Euclid_overall + toc;
    ErrEuclid = 0;
    for a=1:k
        TempClust = find(C1(:,a));
        TempError = N;
        for b=1:k
            TempTempError = length(setdiff(TempClust,TrueClusters{b}));
            TempError = min(TempError,TempTempError);
        end
        ErrEuclid = ErrEuclid + TempError;
    end
    Acc_Euclid_Cluster = 1 - ErrEuclid/N    
   

    % == Update
    Acc_Spectral_Euclid_overall = Acc_Spectral_Euclid_overall + Acc_Euclid_Cluster;

    % ============ 2-weighted Path-distance knn clustering ============ %
    tic
    A2 = ComputeShortestPathAdjMat(Data,15,10,2,sym);
    [C2, ~, ~] = SpectralClustering(A2,k,3);
    Time_Spectral_Path2_overall = Time_Spectral_Path2_overall + toc;

    ErrPath = 0;
    for a=1:k
        TempClust = find(C2(:,a));
        TempError = N;
        for b=1:k
            TempTempError = length(setdiff(TempClust,TrueClusters{b}));
            TempError = min(TempError,TempTempError);
        end
        ErrPath = ErrPath + TempError;
    end
    Acc_Path2_Cluster = 1 - ErrPath/N  

     % == Update vectors
    Acc_Spectral_Path2_overall = Acc_Spectral_Path2_overall + Acc_Path2_Cluster;
    
    % ============ 10-weighted Path-distance knn clustering ============ %
    tic
    A3 = ComputeShortestPathAdjMat(Data,15,10,10,sym);
    [C3, ~, ~] = SpectralClustering(A3,k,3);
    Time_Spectral_LLPD_overall = Time_Spectral_LLPD_overall + toc;

    ErrPath = 0;
    for a=1:k
        TempClust = find(C3(:,a));
        TempError = N;
        for b=1:k
            TempTempError = length(setdiff(TempClust,TrueClusters{b}));
            TempError = min(TempError,TempTempError);
        end
        ErrPath = ErrPath + TempError;
    end
    Acc_Path10_Cluster = 1 - ErrPath/N   

     % == Update vectors
    Acc_Spectral_Path10_overall = Acc_Spectral_Path10_overall + Acc_Path10_Cluster;

     % ============== LLPD knn clustering =========== %
    tic
    A4 = ComputeShortestPathAdjMat_LLPD(Data,15,10,sym);
    Time_Spectral_Path10_overall = Time_Spectral_Path10_overall + toc;
    [C4, ~, ~] = SpectralClustering(A4,k,3);

    ErrLLPD = 0;
    for a=1:k
        TempClust = find(C4(:,a));
        TempError = N;
        for b=1:k
            TempTempError = length(setdiff(TempClust,TrueClusters{b}));
            TempError = min(TempError,TempTempError);
        end
        ErrLLPD = ErrLLPD + TempError;
    end
    Acc_LLPD_Cluster = 1 - ErrLLPD/N    

    % == Update vectors
    Acc_Spectral_LLPD_overall = Acc_Spectral_LLPD_overall + Acc_LLPD_Cluster;
end

Acc_Full_Spectral_Euclid_overall = Acc_Full_Spectral_Euclid_overall./num_trials;
Acc_Spectral_Euclid_overall = Acc_Spectral_Euclid_overall./num_trials;
Acc_Spectral_Path2_overall = Acc_Spectral_Path2_overall./num_trials;
Acc_Spectral_Path10_overall = Acc_Spectral_Path10_overall./num_trials;
Acc_Spectral_LLPD_overall = Acc_Spectral_LLPD_overall./num_trials;

Time_Full_Spectral_Euclid_overall = Time_Full_Spectral_Euclid_overall/num_trials;
Time_Spectral_Path2_overall = Time_Spectral_Path2_overall/num_trials;
Time_Spectral_Path10_overall = Time_Spectral_Path10_overall/num_trials;
Time_Spectral_LLPD_overall = Time_Spectral_LLPD_overall/num_trials;
Time_Spectral_Euclid_overall = Time_Spectral_Euclid_overall/num_trials;




    
    
    



