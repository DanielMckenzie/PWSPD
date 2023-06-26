% ====================================================================== %
% Benchmarking spectral clustering with Euclidean k-nn adjacency matrix vs.
% path distance k-nn adjacency matrix.
% Adjustable to all real data sets
% Daniel Mckenzie
% 1st February 2019
% ======================================================================= %

clear, close all, clc
addpath(genpath('../Utilities'),genpath('../ClusteringAlgorithms'))
addpath(genpath('../Modified_Dijkstra'),genpath('../PathMetric'))
addpath('../../MAT_Files')

% =========================== Parameters ========================= %
colours = ['r', 'b','g','k', 'm', 'c', 'y',];
num_trials = 5;
sym = 'max';

% == COIL
%load('COIL20.mat')
%Points = fea;
%k = 20; n0 = 250; 
%N = 1440;
%y = gnd;

% == USPS
%load('USPS_ready.mat')
%k = 10; n0 = 1100;
%N = k*n0;

% == MNIST
load('MNIST_raw_data.mat')
k = 10;
N = 70000;
Points = images;
y = labels;

% == DrivFace
%load('DrivFace.mat')
%load('DrivFace_Subject.mat')
%k = 4;
%N = 606;
%Points = drivFaceD.data;
%y = subject;

% == OptDigits
%  load('optdigits.mat')
%  Points = double([Xte,Xtr]');
%  y = [Yte,Ytr]';
%  k = 10;
%  N = 5620;

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
    perm = randperm(N);
    [~,permInv] = sort(perm);
    Data = Points(perm,:);
    ynew = y(perm);
    %subjectnew = subject(perm);
    
    % =================== Find the ground truth clusters ============== %
    TrueClusters = cell(k,1);
    % == USPS
    %for a = 1:k
    %    TrueClusters{a} = permInv((a-1)*n0+1:a*n0);
    %end
    
    % == MNIST
   for a = 1:k
       TempClust = find(ynew == a-1);
       TrueClusters{a} = TempClust;
   end
   disp('Loaded the MNIST data!')
   
    % == DrivFace
    %for a = 1:k
    %   TempClust = find(ynew == a);
    %   TrueClusters{a} = TempClust;
    %end
    
    % == OptDigits
%     for a = 1:k
%         TempClust = find(ynew == a-1);
%         TrueClusters{a} = TempClust;
%     end
    
    % == COIL20
%     for a = 1:k
%         TempClust = find(ynew == a);
%         TrueClusters{a} = TempClust;
%     end
    
     % ============== Full matrix Euclidean clustering =========== %
%     tic
%     A0 = SimilarityMatrixZMP(Data,10);
%     [C0, ~, ~] = SpectralClustering(A0, k, 3);
%     Time_Full_Spectral_Euclid_overall = Time_Full_Spectral_Euclid_overall + toc;
%     ErrFullEuclid = 0;
%     for a=1:k
%         TempClust = find(C0(:,a));
%         TempError = N;
%         for b=1:k
%             TempTempError = length(setdiff(TempClust,TrueClusters{b}));
%             TempError = min(TempError,TempTempError);
%         end
%         ErrFullEuclid = ErrFullEuclid + TempError;
%     end
%     Acc_Full_Euclid_Cluster = 1 - ErrFullEuclid/N 
%     
%      % == Update
%     Acc_Full_Spectral_Euclid_overall = Acc_Full_Spectral_Euclid_overall + Acc_Full_Euclid_Cluster;
   
     
    % ============== Euclidean knn clustering =========== %
    tic
    A1 = CreateKNN_Max_from_Data(Data,15,10);
    disp('Created A1!')
    [C1, ~, ~] = SpectralClustering(A1, k, 3);
    disp('Finished Spectral Clustering with A1')
    Time_Spectral_Euclid_overall = Time_Spectral_Euclid_overall + toc
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
    clear A1 C1

    % ============ 2-weighted Path-distance knn clustering ============ %
    tic
    A2 = ComputeShortestPathAdjMat(Data,15,10,2,sym);
    disp('Created A2!')
    [C2, ~, ~] = SpectralClustering(A2,k,3);
    disp('Finished Spectral Clustering with A2')
    Time_Spectral_Path2_overall = Time_Spectral_Path2_overall + toc

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
    clear A2 C2
    
    % ============ 10-weighted Path-distance knn clustering ============ %
    tic
    A3 = ComputeShortestPathAdjMat(Data,15,10,10,sym);
    disp('created A3!')
    [C3, ~, ~] = SpectralClustering(A3,k,3);
    disp('Finished Spectral Clustering with A3')
    Time_Spectral_Path10_overall = Time_Spectral_Path10_overall + toc;

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
    clear A3 C3

     % ============== LLPD knn clustering =========== %
    tic
    A4 = ComputeShortestPathAdjMat_LLPD(Data,15,10,sym);
    disp('created A4')
    [C4, ~, ~] = SpectralClustering(A4,k,3);
    disp('Finished Spectral Clustering with A4')
    Time_Spectral_LLPD_overall = Time_Spectral_LLPD_overall + toc;
    
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
    clear A4 C4
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






    
    
    



