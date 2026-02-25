tic
clear all
format long

%%
addpath(genpath('../GenLouvain-master')); %Use the GenLouvain directory to the path

load('data/groupwise_oscillating_neurons_I=5_combined_comm=6_2_new.mat')
groups = load('data/sbm_partition_I=5_combined_comm=6_2_new.mat');
groupIndices = groups.groupIndices;
% load('data/signaling_group_wise_oscillating_neurons_I=5_combined.mat')
% groups = load('data/signalingCommunities.mat');
% groupIndices = groups.signalling_comm_indices;

% % groupIndices is a 6x1 cell, each containing a vector
% sizes = cellfun(@numel, groupIndices);   % get sizes of each group
% 
% [~, sortIdx] = sort(sizes, 'descend');   % sort by size (largest first)
% % [~, sortIdx] = sort(sizes, 'ascend');  % use this for smallest first
% 
% groupIndices = groupIndices(sortIdx);

nGroups = numel(groupIndices);

I_set = [50000;75000;100000;200000;300000]; %External current
perturbation = 2.5e-3; % Perturbation from the equilibrium
end_time = 200; % End time

threshold_set = [0.55 0.6 0.65 0.7];%0.4:0.05:0.9;  % set your threshold for Louvain Algorithm

Plotting.bool = false;
%%
state = cell(6,5,length(threshold_set));
grp_wise_avg_kuramoto_order = cell(6,1);
for i=1:6
    grp_wise_avg_kuramoto_order{i} = zeros(6,6);
end
count = zeros(6,1);
for m = 1:5
    I = I_set(m);
    B_I = B(:,m);
    for k =1:nGroups
        dummy = B_I{k};
        for l = 1:length(dummy)
            fprintf('I = %d, group = %d, neuron = %d \n', I, k, l)
            neuron = dummy(l); %Stimulated neuron
            Results = Run_Model(neuron,I,end_time,Plotting,cd,perturbation);
    
            t = Results.time;
            idx = 1;%length(t)/10;
            V = Results.V(:,idx:end)';
            S = Results.S(:,idx:end)';
            t = Results.time(idx:end);
            
            % Center and normalize signals
            s_norm = Centered_Signal(S);
            v_norm = Centered_Signal(V);
            
            % phi = atan2(S,V);
            phi = atan2(s_norm,v_norm); % Geometric phase.
            phi(isnan(phi)) = 0;
            
            % Preallocate matrix for average Kuramoto order parameter
            Rmat = zeros(nGroups);
            
            for i = 1:nGroups
                for j = 1:nGroups
                    % neurons in group i and j
                    neurons = unique([groupIndices{i}; groupIndices{j}]);
                    
                    % extract their phases (nTime x numNeurons)
                    theta = phi(:, neurons);
                    
                    % compute Kuramoto order parameter at each time
                    R_t = abs(mean(exp(1i*theta), 2));
                    
                    % average over time
                    Rmat(i,j) = mean(R_t);
                end
            end
            grp_wise_avg_kuramoto_order{k} = grp_wise_avg_kuramoto_order{k} + Rmat;
            count(k) = count(k) + 1;
            for n = 1:length(threshold_set)
                threshold = threshold_set(n);
                % A_filtered = A(P(:,1),P(:,1));
                R_bin = double(Rmat >= threshold);
                
                [S,Q] = genlouvain(R_bin); %Use the Louvain.m here for community assignment
                state{k,m,n} = [state{k,m,n};S'];
                % if isscalar(unique(S))
                %     state(k,1,m,n) = state(k,1,m,n) + 1;
                % elseif length(unique(S)) == nGroups
                %     state(k,2,m,n) = state(k,2,m,n) + 1;
                % else
                %     state(k,3,m,n) = state(k,3,m,n) + 1;
                % end
            end
        end
    end
end
for i = 1:6
    grp_wise_avg_kuramoto_order{i} = grp_wise_avg_kuramoto_order{i}./count(i);
end
save('data/state_patterns_th=0.55_0.6_0.65_0.7_I=5_combined_comm=6_2_new.mat', 'state');
save('data/driving_group_wise_avg_kuramoto_order_comm=6_2_new.mat', 'grp_wise_avg_kuramoto_order');

%%
toc






%%
%%%%%%%%%%%%%callable functions used in the above code****************
function Signal_NORM = Centered_Signal(Signal)
    Signal_MAX = max(Signal);
    Signal_MIN = min(Signal);
    Signal_AMPLITUDE = (Signal_MAX - Signal_MIN)/2;
    Signal_SHIFT = Signal_MAX - Signal_AMPLITUDE;
    Signal_TEMP = Signal - Signal_SHIFT;
    Signal_NORM = Signal_TEMP./Signal_AMPLITUDE;
end


function x_analytic = custom_hilbert(x)
    % Compute the Hilbert transform (analytic signal) of input signal x
    % Equivalent to MATLAB's hilbert(x) from the Signal Processing Toolbox
    %
    % Input:
    %   x : Real-valued signal (vector or matrix)
    % Output:
    %   x_analytic : Analytic signal (complex signal = x + j*HilbertTransform(x))

    % Get the size of the input
    N = size(x, 1);

    % Compute the FFT of the input signal along the first dimension
    X = fft(x, [], 1);

    % Create the phase shift vector for the analytic signal
    h = zeros(size(X));
    if mod(N, 2) == 0
        % Even length
        h(1) = 1;
        h(2:N/2) = 2;
        h(N/2 + 1) = 1; % Nyquist frequency
    else
        % Odd length
        h(1) = 1;
        h(2:(N + 1)/2) = 2;
    end

    % Apply the phase shift to the FFT
    X_analytic = X .* h;

    % Compute the inverse FFT to get the analytic signal
    x_analytic = ifft(X_analytic, [], 1, 'symmetric');
end

function neuron_number = label2num(neuron_label,type_neuron_label)
    neuron_number = zeros(length(type_neuron_label),1);
    for i = 1:length(type_neuron_label)
        neuron_number(i) = find(neuron_label == type_neuron_label(i));
    end
end