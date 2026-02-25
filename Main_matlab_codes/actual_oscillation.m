tic
clear all
format long
data = importdata("Worm279dir.mat"); %load the neuron data
Neuron_Labels = data.Worm279_labelled;
N = 279; % Total number of neurons
unstable_neurons = []; % Create a blank array for unstable neurons

I = 275000; % Change the external current here
perturbation = 2.5e-3; % Perturbation from the equilibrium
end_time = 200; % End time

% This loop finds the neurons that are unstable to the external current
for i=1:N
    fprintf('Finding unstable neurons: %.2f%%\n', i/N*100)
    EV = zeros(2*N,1);
    EV = Find_Hopf_Bifurcation_Stepper(I,i,cd);
    if any(real(EV(:)) > 0)
        unstable_neurons = [unstable_neurons; i];
    end
end

% Save the unstable neuron data with the neuron number, name, and category
labels = Neuron_Labels(unstable_neurons,1);
type = Neuron_Labels(unstable_neurons,2);
unstable_neurons_new = [num2cell(unstable_neurons(:)), labels(:), type(:)];
save('data/unstable-neurons-for-I=275000-tr-weighted.mat', 'unstable_neurons_new');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now look for the oscillating neurons
x = unstable_neurons_new;
x = cell2mat(x(:,1));
oscillation = zeros(length(x),3); % store mean oscillation

% This loop calculates the mean oscillation of the unstable neurons
for i = 1:length(x)
    fprintf('Finding actual oscillation: %.2f%%\n', i/length(x)*100)
    neuron = x(i); %Stimulated neuron
    Plotting.bool = false;
    Results = Run_Model(neuron,I,end_time,Plotting,cd,perturbation);
    
    t = Results.time;
    idx = length(t)/10;
    V = Results.V(:,idx:end)';
    S = Results.S(:,idx:end)';
    t = Results.time(idx:end);
    
    a_v=max(V);
    b_v=min(V);
    c_v=a_v-b_v;
    a_s=max(S);
    b_s=min(S);
    c_s=a_s-b_s;
    oscillation(i,1) = neuron;
    oscillation(i,2) = mean(c_v);
    oscillation(i,3) = mean(c_s);
end
oscillating_neurons_new = oscillation(round(oscillation(:,2),4) ~= 0, :); % filter out the oscillating neuron
save('data/oscillating-neurons-for-I=275000-tr-weighted.mat', 'oscillating_neurons_new'); % save the data
toc