clear all
clc
addpath(genpath('../analysis tools'));
% --- inputs ---
groups = load('data/sbm_partition_I=5_combined_comm=6_3.mat');
% groups = load('data/sbm_partition_I=11_combined_comm=6.mat');
groupIndices = groups.groupIndices;
refPart  = groupIndices;   % k×1 cell
numRuns  = 50;

k = numel(refPart);

% assume node IDs are integers from 1..N
N = max(cellfun(@max, refPart));   % adjust if needed

% nodeGroups(node, run) = group index (1..k) after alignment
nodeGroups = zeros(N, numRuns);

% perms(i, run) = which group in the raw partition matched reference group i
perms = zeros(k, numRuns);

A1 = load('data/co_occurrence_matrix_I=5_combined_comm=8.mat');
A = A1.G;
A_org = A;
A = A - mean(A(:)); A = A./max(A(:));

A(1:(N + 1):end) = 0;               % remove diagonal

W_distr = 'normal';                 % distribution over edge weights
E_distr = 'DC';                     % Degree Corrected
alph = 0.5;                        % balance weight/edge probability distributions (alph = 0 means we conly care about weights)
muMaxIter = 250;                    % number of optimization steps for mu parameter
mainMaxIter = 250;                  % number of optimization steps for main loop
mainTol = 1.0000e-03;               % tolerance -- smaller = better fit, but longer runtime
muTol = 1.0000e-03;                 % tolerance -- smaller = better fit, but longer runtime
edgeList = Adj2Edg(A);              % convert matrix to edges

for r = 1:numRuns
    r
    % ========================
    % 1. Run your partition method
    % ========================
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
    part = groupIndices;

    % ========================
    % 2. Build similarity matrix with reference
    %    S(iRef, jRun) = overlap size between ref group iRef and group jRun
    % ========================
    S = zeros(k);
    for iRef = 1:k
        refNodes = refPart{iRef};
        for jRun = 1:k
            runNodes = part{jRun};
            % similarity = size of intersection
            S(iRef, jRun) = numel(intersect(refNodes, runNodes));

            % If you prefer Jaccard similarity, use this instead:
            % S(iRef, jRun) = numel(intersect(refNodes, runNodes)) / ...
            %                numel(union(refNodes, runNodes));
        end
    end

    % ========================
    % 3. Find best matching between reference groups and run groups
    % ========================
    perm = greedy_best_matching(S);   % perm(iRef) = index in 'part'
    perms(:, r) = perm(:);

    % ========================
    % 4. Reorder the groups of this run to match the reference labels
    % ========================
    alignedPart = part(perm);   % alignedPart{iRef} ~ refPart{iRef}

    % ========================
    % 5. Record, for each node, which (aligned) group it belongs to
    % ========================
    for g = 1:k
        nodesInG = alignedPart{g};
        nodeGroups(nodesInG, r) = g;
    end
end

% ---------------------------------------------------------
% After the loop:
% - nodeGroups(n,:) is the sequence of partition labels (1..k)
%   that node n belongs to across the 10 runs (after alignment).
% - unique(nodeGroups(n,:)) is the *set* of partitions it belongs to.
% ---------------------------------------------------------
function perm = greedy_best_matching(S)
% S: k×k similarity matrix, rows = reference groups, columns = run groups
    k = size(S, 1);
    perm = zeros(1, k);
    usedCols = false(1, k);

    for i = 1:k
        scores = S(i, :);
        scores(usedCols) = -inf;   % don’t reuse groups
        [~, jBest] = max(scores);
        perm(i) = jBest;
        usedCols(jBest) = true;
    end
end
%%
[N, numRuns] = size(nodeGroups);
k = max(nodeGroups(:));   % number of groups

mode1 = zeros(N,1);   % first mode (most frequent label)
mode2 = zeros(N,1);   % second mode (2nd most frequent label)

for n = 1:N
    labels = nodeGroups(n,:);   % labels of node n across runs

    % Unique labels and their counts
    u = unique(labels);
    counts = zeros(size(u));

    for idx = 1:numel(u)
        counts(idx) = sum(labels == u(idx));
    end

    % Sort by frequency (descending)
    [~, sortIdx] = sort(counts, 'descend');
    sortedLabels = u(sortIdx);

    % First mode: always exists
    mode1(n) = sortedLabels(1);

    % Second mode: if there is a second distinct label
    if numel(sortedLabels) >= 2
        mode2(n) = sortedLabels(2);
    else
        % If node only ever appears in one group, there is no 2nd mode
        % Option 1 (default here): mark as NaN and ignore later in 2nd partition
        mode2(n) = mode1(n);

        % Option 2 (if you want every node in both partitions),
        % you can instead do:
        % mode2(n) = sortedLabels(1);
    end
end

% Initialize partitions
part_mode1 = cell(k,1);
part_mode2 = cell(k,1);

% Partition from first modes
for g = 1:k
    part_mode1{g} = find(mode1 == g);     % nodes whose primary label is g
end

% Partition from second modes
for g = 1:k
    part_mode2{g} = find(mode2 == g);     % nodes whose secondary label is g
end

save('data/sbm_partition_I=5_combined_comm=6_50runs_two_modes_alpha_0.5.mat', 'part_mode1', 'part_mode2');