function dy = Conductance_ODE(y,N,I,Parameters)
    % The Neuronal Model being used. A conductance model.
    % Arguments:
    % ----------
    % y: Nx2 vector; 1:N Voltage, N+1:2*N Synaptic Activity Variable
    % N: Number of Neurons, should be N = 279 for the C. Elegans.
    % I: Nx1 vector; The applied current.
    % Parameters: Struct; All the parameters the model uses.
    
    %%% Parameters:
    tau = Parameters.tau;     %  s; Synaptic Time Delay Constant.
    G_L = Parameters.G_L;     % pS; Membrane Conductance.
    g = Parameters.g;         % pS; Synapse/Gap Conductance.
    V_L = Parameters.V_L;     % mV; Leakage Potential.
    L_gap = Parameters.L_gap; % Gap junction adjacency Matrix.
    W_syn = Parameters.W_syn; % Synaptic Adjacency Matrix.
    E = Parameters.E;         % Reversal Potentials.
    a_r = Parameters.a_r;     %  s; Time Rise Constant.
    a_d = Parameters.a_d;     %  s; Time Delay Constant.
    beta = Parameters.beta;   % mV; Sigmoid Width.
    V_th = Parameters.V_th;   % mV; Threshold Voltage.

    
    %%% Variables
    V = y(1:N,1);     % mV; Voltage
    S = y(N+1:2*N,1); % Synaptic Activity Variable.

    %%% Apply Noise
    noise_var = 0.00005; 
    noise_V = noise_var * randn(N,1);
    noise_S = noise_var * randn(N,1);
    
    %%% ODE
    V_mem = -(G_L/g).*(V-V_L);
    V_gap = L_gap*V;
    V_syn = V.*(W_syn'*S) - W_syn'*(S.*E);
    dy(1:N,1) = (V_mem - V_gap - V_syn + I)/tau + noise_V;
%     dy(1:N,1) = (V_mem - V_gap - V_syn)/tau + I + noise_V;
    
    phi = 1./(1+exp(-beta*(V - V_th)));
    dy(N+1:2*N,1) = a_r*(phi.*(1-S)) - a_d*S + noise_S;
    
end