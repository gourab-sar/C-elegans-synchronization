clear all
% close all
% clc
% addpath(genpath('../WSBM_v1.2'));      % add wsbm toolbox to path (http://tuvalu.santafe.edu/~aaronc/wsbm/)
% Add Analysis Tool Dir
addpath(genpath('../analysis tools'));

% A1 = load('data/co_occurrence_matrix_I=5_combined_comm=8.mat');
A1 = load('data/co_occurrence_matrix_I=5_combined_comm=4_alpha_05.mat');
A = A1.G;
A_org = A;
A = A - mean(A(:)); A = A./max(A(:));

% A1 = load('data/pearson_correlation-matrix-for-I=50000-tr-weighted.mat');
% A1 = A1.v_pearson_corr_matrix;
% A = A1(:,:,23);
% A_org = A;

N = length(A);                      % number of nodes
A(1:(N + 1):end) = 0;               % remove diagonal

W_distr = 'normal';                 % distribution over edge weights
E_distr = 'DC';                     % Degree Corrected
alph = 0.5;                        % balance weight/edge probability distributions (alph = 0 means we conly care about weights)
muMaxIter = 250;                    % number of optimization steps for mu parameter
mainMaxIter = 250;                  % number of optimization steps for main loop
mainTol = 1.0000e-03;               % tolerance -- smaller = better fit, but longer runtime
muTol = 1.0000e-03;                 % tolerance -- smaller = better fit, but longer runtime

k = 4;                              % number of communities to detect
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
        'mainTol',mainTol,...
        'muTol',muTol,...
        'alpha',alph);
    mx = max(length(unique(ci)));   % number of unique communities
    
    toc
end
[C,Morph] = fcn_comm_motifs(A,ci);  % get community motif types -- type 'help fcn_comm_motifs' in command line for more information

[~, sortIdx] = sort(ci);
group_sorted = ci(sortIdx);
A_grouped = A_org(sortIdx, sortIdx);

group = ci;
% group = 100x1 vector of group indices
uniqueGroups = unique(group);   % distinct groups
groupLengths = zeros(length(uniqueGroups),1);
groupIndices = cell(length(uniqueGroups),1);
for kk = 1:length(uniqueGroups)
    g = uniqueGroups(kk);
    groupLengths(kk) = length(find(group == g));
    groupIndices{kk} = find(group == g);   % indices of neurons in group g
end
mingroupLength = min(groupLengths);


B1 = groupIndices;
orderedIdx = vertcat(B1{:});   % concatenate indices of all groups
% 3) Compute group sizes and boundaries
counts = cellfun(@numel, B1);  % number of neurons per group
boundaries = cumsum(counts);

% 4) Plot with boundaries
figure;
% Adjust figure window to stretch horizontally
screenSize = get(0, 'ScreenSize'); % Get screen dimensions
figureWidth = screenSize(3) * 0.8; % Full screen width
figureHeight = screenSize(4) * 0.4; % 60% of screen height (adjust as needed)
figureLeft = (screenSize(3) - figureWidth)/2; % Start at left edge
figureBottom = (screenSize(4) - figureHeight)/2; % Center vertically (optional)
set(gcf, 'Position', [figureLeft figureBottom figureWidth figureHeight]);

subplot(1,2,1)
imagesc(A_org)
subplot(1,2,2)
imagesc(A_grouped);
colormap(parula);
colorbar;
axis square;

hold on;
for b = boundaries(1:end-1)   % skip last (matrix edge)
    xline(b + 0.5, 'r', 'LineWidth', 2);
    yline(b + 0.5, 'r', 'LineWidth', 2);
end
hold off;

% save('data/sbm_partition_I=5_combined_comm=4.mat', 'groupIndices');