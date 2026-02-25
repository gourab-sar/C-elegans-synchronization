clear all

load('data/groupwise_oscillating_neurons_I=5_combined_comm=6_2_new.mat');
% load('data/signaling_group_wise_oscillating_neurons_I=5_combined.mat');

B_union_ordered = cell(6,1);
for i = 1:6
    for j = 1:5
        B_union_ordered{i} = [B_union_ordered{i};B{i,j}];
    end
end

load('data/state_patterns_th=0.55_0.6_0.65_0.7_I=5_combined_comm=6_2_new.mat');
% load('data/state_patterns_signaling_groups_th=0.55_0.6_0.65_0.7_I=5_combined.mat');
state_55 = state(:,:,1);
state_60 = state(:,:,2);
state_65 = state(:,:,3);
state_70 = state(:,:,4);

state_60_union = cell(6,1);
for i = 1:6
    for j = 1:5
        state_60_union{i} = [state_60_union{i};state_60{i,j}];
    end
end

state_65_union = cell(6,1);
for i = 1:6
    for j = 1:5
        state_65_union{i} = [state_65_union{i};state_65{i,j}];
    end
end

% Get unique rows and their first occurrence indices
[UniqueCombinations, ~, ic] = unique(state_65_union{4}, 'rows', 'stable');

% Count occurrences
Counts = accumarray(ic, 1);

% Show results
table(UniqueCombinations, Counts)


numUnique = size(UniqueCombinations, 1);

%% combine all neurons to their patterns
assocCells = cell(6,1); % output, same structure

for g = 1:6
    neurons = B_union_ordered{g};      % column/row vector of indices
    vectors = state_65_union{g};      % matrix (#neurons x 6)

    % Ensure sizes match
    if numel(neurons) ~= size(vectors,1)
        error('Group %d: neuron count does not match number of vectors!', g);
    end

    % Store as Nx2 cell array: {neuronIndex, correspondingVector}
    assocCells{g} = [num2cell(neurons(:)) num2cell(vectors,2)];
end
%% find unique patterns
uniqueAssoc = cell(6,1);

for g = 1:6
    pairs = assocCells{g};
    neurons = cell2mat(pairs(:,1));
    vectors = cell2mat(pairs(:,2));

    % Combine neuron index and vector into one row
    combined = [neurons(:) vectors];
    
    % Take unique rows
    [~, ia] = unique(combined, 'rows');
    
    % Rebuild unique assocCells
    uniqueAssoc{g} = [num2cell(combined(ia,1)) num2cell(combined(ia,2:end),2)];
end

%%
state_unique = cell(4,1);
for k =1:4
    % combine all neurons to their patterns
    assocCells = cell(6,1); % output, same structure

    state_union = cell(6,1);
    for i = 1:6
        for j = 1:5
            state_union{i} = [state_union{i};state{i,j,k}];
        end
    end
    
    for g = 1:6
        neurons = B_union_ordered{g};      % column/row vector of indices
        vectors = state_union{g};      % matrix (#neurons x 6)
    
        % Ensure sizes match
        if numel(neurons) ~= size(vectors,1)
            error('Group %d: neuron count does not match number of vectors!', g);
        end
    
        % Store as Nx2 cell array: {neuronIndex, correspondingVector}
        assocCells{g} = [num2cell(neurons(:)) num2cell(vectors,2)];
    end
    
    % find unique patterns
    uniqueAssoc = cell(6,1);
    
    for g = 1:6
        pairs = assocCells{g};
        neurons = cell2mat(pairs(:,1));
        vectors = cell2mat(pairs(:,2));
    
        % Combine neuron index and vector into one row
        combined = [neurons(:) vectors];
        
        % Take unique rows
        [~, ia] = unique(combined, 'rows');
        
        % Rebuild unique assocCells
        uniqueAssoc{g} = [num2cell(combined(ia,1)) num2cell(combined(ia,2:end),2)];
    end
    state_unique{k} = uniqueAssoc;
end

save('data/state_patterns_th=0.55_0.6_0.65_0.7_I=5_combined_comm=6_2_new.mat', 'state_unique', '-append');
% save('data/state_patterns_signaling_groups_th=0.55_0.6_0.65_0.7_I=5_combined.mat', 'state_unique', '-append');


%%
num_thresh = numel(state_unique);
num_groups = numel(state_unique{1});

all_vectors = cell(num_thresh, num_groups);

for t = 1:num_thresh
    for g = 1:num_groups
        group_data = state_unique{t}{g};   % m×2 cell
        vectors = group_data(:,2);         % column 2 = vectors
        all_vectors{t,g} = cell2mat(vectors); 
        % cell2mat works if all vectors have same length
        % otherwise, keep as cell array
    end
end
state_unique_python = all_vectors;
save('data/state_patterns_th=0.55_0.6_0.65_0.7_I=5_combined_comm=6_2_new.mat', 'state_unique_python', '-append');
% save('data/state_patterns_signaling_groups_th=0.55_0.6_0.65_0.7_I=5_combined.mat', 'state_unique_python', '-append');



% state_65_all = [];
% for i=1:6
%     state_65_all = [state_65_all;state_65_union{i}];
% end
% 
% % Get unique rows and their first occurrence indices
% [UniqueCombinations, ~, ic] = unique(state_65_all, 'rows', 'stable');
% 
% % Count occurrences
% Counts = accumarray(ic, 1);
% 
% % Show results
% table(UniqueCombinations, Counts)