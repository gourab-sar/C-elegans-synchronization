clear all
close all
clc
% addpath(genpath('../WSBM_v1.2'));      % add wsbm toolbox to path (http://tuvalu.santafe.edu/~aaronc/wsbm/)
% Add Analysis Tool Dir
addpath(genpath('../analysis tools')); 

% Install MEX Files (Optional)
%InstallMEXFiles


A1 = load('data/pearson_correlation-matrix-for-I=50000-tr-weighted.mat');
A1 = A1.v_pearson_corr_matrix;
A2 = load('data/pearson_correlation-matrix-for-I=75000-tr-weighted.mat');
A2 = A2.v_pearson_corr_matrix;
A3 = load('data/pearson_correlation-matrix-for-I=100000-tr-weighted.mat');
A3 = A3.v_pearson_corr_matrix;
A4 = load('data/pearson_correlation-matrix-for-I=125000-tr-weighted.mat');
A4 = A4.v_pearson_corr_matrix;
A5 = load('data/pearson_correlation-matrix-for-I=150000-tr-weighted.mat');
A5 = A5.v_pearson_corr_matrix;
A6 = load('data/pearson_correlation-matrix-for-I=175000-tr-weighted.mat');
A6 = A6.v_pearson_corr_matrix;
A7 = load('data/pearson_correlation-matrix-for-I=200000-tr-weighted.mat');
A7 = A7.v_pearson_corr_matrix;
A8 = load('data/pearson_correlation-matrix-for-I=225000-tr-weighted.mat');
A8 = A8.v_pearson_corr_matrix;
A9 = load('data/pearson_correlation-matrix-for-I=250000-tr-weighted.mat');
A9 = A9.v_pearson_corr_matrix;
A10 = load('data/pearson_correlation-matrix-for-I=275000-tr-weighted.mat');
A10 = A10.v_pearson_corr_matrix;
A11 = load('data/pearson_correlation-matrix-for-I=300000-tr-weighted.mat');
A11 = A11.v_pearson_corr_matrix;

% AA = cat(3, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11);
AA = cat(3, A1, A2, A3, A7, A11);

%%
N = size(AA,1);   %number of neurons

% group = 279x1 vector of group indices
G = zeros(N);

W_distr = 'normal';                 % distribution over edge weights
E_distr = 'DC';                     % Degree Corrected
alph = 0.5;                           % balance weight/edge probability distributions (alph = 0 means we conly care about weights)
muMaxIter = 250;                    % number of optimization steps for mu parameter
mainMaxIter = 250;                  % number of optimization steps for main loop
mainTol = 1.0000e-03;               % tolerance -- smaller = better fit, but longer runtime
muTol = 1.0000e-03;                 % tolerance -- smaller = better fit, but longer runtime

k = 4;                              % number of communities to detect

for l = 1:size(AA,3)
    l
    A= AA(:,:,l);
    A(1:(N + 1):end) = 0;               % remove diagonal
    
    edgeList = Adj2Edg(A);              % convert matrix to edges
    mx = inf;                           % initialize number of communities to inf
    mingroupLength = 0;
    while mx~=k %~(mx == k && mingroupLength >= 4)                      % run algorithm until you get k communities -- due to "tie-breaking" you can get fewer than k
        
        tic
        [ci,m] = wsbm(edgeList,k,...    % run model
            'NumTrials',1,...
            'muMaxIter',muMaxIter,...
            'mainMaxIter',mainMaxIter,...
            'W_distr',W_distr,...
            'E_distr',E_distr,...
            'mainTol',mainTol,...
            'muTol',muTol,...
            'alpha',alph);
        mx = max(length(unique(ci)));   % number of unique communities
        % group = ci;
        % uniqueGroups = unique(group);   % distinct groups
        % groupLengths = zeros(length(uniqueGroups),1);
        % for kk = 1:length(uniqueGroups)
        %     g = uniqueGroups(kk);
        %     groupLengths(kk) = length(find(group == g));   % indices of neurons in group g
        % end
        % mingroupLength = min(groupLengths);
        toc
    end
    [C,Morph] = fcn_comm_motifs(A,ci);  % get community motif types -- type 'help fcn_comm_motifs' in command line for more information
    
    [~, sortIdx] = sort(ci);
    group_sorted = ci(sortIdx);
    A_grouped = A(sortIdx, sortIdx);
    
    group = ci;
    
    
    for i = 1:N
        for j = 1:N
            if group(i) == group(j)
                G(i,j) = G(i,j)+1;
            end
        end
    end
end
save('data/co_occurrence_matrix_I=5_combined_comm=4_alpha_05.mat', 'G');