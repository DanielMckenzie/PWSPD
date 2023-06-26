% ====================================================================== %
% Testing power weighted path clustering for various powers.
% All three synthetic data sets.
% Daniel Mckenzie
% 25th January 2019
% ====================================================================== %

clear, close all, clc
addpath(genpath('../Utilities'),genpath('../ClusteringAlgorithms'))
addpath(genpath('../Modified_Dijkstra'),genpath('../PathMetric'))
addpath('../MAT_Files')

% =========================== Parameters ========================= %
colours = ['r', 'b','g','k', 'm', 'c', 'y',];
n0 = 300; noise_level = 0.15; k=3;
N = k*n0;

num_p = 20;
num_trials = 10;

% == Lines
%ysep=1; len = 5;

% == Moons
r1 = 1;
r2 = 1;
r3 = 1.5;



% ==================== Define the array of interest ================ %
Acc_Spectral_path_vec = zeros(num_p,3);

for m = 1:3 % Iterate over the three ambient dims to try
    if m == 1
        ambient_dim = 10;
    elseif m ==2
        ambient_dim = 50;
    else
        ambient_dim = 100;
    end
    
    for j = 1:num_trials
         % ==================== Generate and Permute the data ============== %
        % == Lines
        %[Points2D,Points] = Generate3Lines(ysep,len,n0,noise_level,ambient_dim);
         % == Moons
        [Points2D,Points] = Generate3Moons(r1,r2,r3,n0,noise_level,ambient_dim);
   
        perm = randperm(N);
        [~,permInv] = sort(perm);
        Data = Points(perm,:);

        % =================== Find the ground truth clusters ============== %
        TrueClusters = cell(k,1);
        for a = 1:k
            TrueClusters{a} = permInv((a-1)*n0+1:a*n0);
        end

        for i = 1:num_p
            % ========= Build appropriate power-weighted adj mat =========== %
            A2 = ComputeShortestPathAdjMat(Data,15,10,i,'max');
            try 
                [C2, ~, ~] = SpectralClustering(A2,k,3);
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
                Acc_Spectral_path_vec(i,m) = Acc_Spectral_path_vec(i,m) + 1 - ErrPath/N
            catch
                disp(['Error Caught! p = ',num2str(i)])
                Acc_Spectral_path_vec(i,m) = Acc_Spectral_path_vec(i,m);
            end
        end
    end
end

Acc_Spectral_path_vec = 100*Acc_Spectral_path_vec./num_trials;

% ================ Plot resulting data ================== %
plot(1:num_p, Acc_Spectral_path_vec(:,1),'LineWidth',2);
hold on
plot(1:num_p, Acc_Spectral_path_vec(:,2),'LineWidth',2);
plot(1:num_p, Acc_Spectral_path_vec(:,3),'LineWidth',2);
legend({'D=10','D=50','D=100'},'FontSize',18)
set(gca,'FontSize',14)


    