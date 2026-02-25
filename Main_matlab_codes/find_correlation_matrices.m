tic
clear all
format long

N =279; % Total number of neurons
I = 275000; %External current
perturbation = 2.5e-3; % Perturbation from the equilibrium
end_time = 200; % End time

osc = load('data/oscillating-neurons-for-I=275000-tr-weighted.mat');
osc = osc.oscillating_neurons_new; osc = osc(:,1);
n_osc = length(osc);

v_pearson_corr_matrix = zeros(N,N,n_osc);
v_norm_pearson_corr_matrix = zeros(N,N,n_osc);
s_pearson_corr_matrix = zeros(N,N,n_osc);
s_norm_pearson_corr_matrix = zeros(N,N,n_osc);
phi_pearson_corr_matrix = zeros(N,N,n_osc);

for i = 1:n_osc
    fprintf('Calculating correlation matrix: %d of %d %.2f%%\n', i, n_osc, i/n_osc*100)
    neuron = osc(i,1); %Stimulated neuron
    Plotting.bool = false;
    Results = Run_Model(neuron,I,end_time,Plotting,cd,perturbation);
    
    t = Results.time;
    idx = length(t)/10;
    V = Results.V(:,idx:end)';
    S = Results.S(:,idx:end)';
    t = Results.time(idx:end);
    
    % Center and normalize signals
    s_norm = Centered_Signal(S);
    v_norm = Centered_Signal(V);
    
    phi = atan2(s_norm,v_norm); % Geometric phase.
%     phi(isnan(phi)) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%CORRELATION%%%%%%%%%%%%
    keep_idx = setdiff(1:1:279,[122,123,149,154]);
    V1 = V;%(:,keep_idx);
    S1 = S;%(:,keep_idx);
    v_norm1 = v_norm;%(:,keep_idx);
    s_norm1 = s_norm;%(:,keep_idx);
    phi1 = phi;%(:,keep_idx);

    v_org_corr = corr(V1);
    v_org_corr(isnan(v_org_corr)) = 0;
    v_org_corr = round(v_org_corr,4);

    v_norm_corr = corr(v_norm1);
    v_norm_corr(isnan(v_norm_corr)) = 0;
    v_norm_corr = round(v_norm_corr,4);

    s_org_corr = corr(S1);
    s_org_corr(isnan(s_org_corr)) = 0;
    s_org_corr = round(s_org_corr,4);

    s_norm_corr = corr(s_norm1);
    s_norm_corr(isnan(s_norm_corr)) = 0;
    s_norm_corr = round(s_norm_corr,4);

    phi_corr = corr(phi1);
    phi_corr(isnan(phi_corr)) = 0;
    phi_corr = round(phi_corr,4);

    v_pearson_corr_matrix(:,:,i) = v_org_corr;
    v_norm_pearson_corr_matrix(:,:,i) = v_norm_corr;
    s_pearson_corr_matrix(:,:,i) = s_org_corr;
    s_norm_pearson_corr_matrix(:,:,i) = s_norm_corr;
    phi_pearson_corr_matrix(:,:,i) = phi_corr;

end
 
save('data/pearson_correlation-matrix-for-I=275000-tr-weighted.mat', 'v_pearson_corr_matrix',...
    'v_norm_pearson_corr_matrix', 's_pearson_corr_matrix', 's_norm_pearson_corr_matrix',...
    'phi_pearson_corr_matrix');


toc

function Signal_NORM = Centered_Signal(Signal)
    Signal_MAX = max(Signal);
    Signal_MIN = min(Signal);
    Signal_AMPLITUDE = (Signal_MAX - Signal_MIN)/2;
    Signal_SHIFT = Signal_MAX - Signal_AMPLITUDE;
    Signal_TEMP = Signal - Signal_SHIFT;
    Signal_NORM = Signal_TEMP./Signal_AMPLITUDE;
end
