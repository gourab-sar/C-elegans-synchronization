%%% Running the Conductance Model
function Results = Run_Model(Perturbed_Neurons,I_amplitude,end_time,Plotting,path,perturb)
    %% Parameters:
    N = 279; % Total Number of Neurons.
%     Network_Parameters = importdata(path+"/Worm279dir.mat");
    Network_Parameters = importdata("C_Elegans_Network_Data.mat");
    Reversal_Potential = importdata(path+"/Reversal_Potential.mat");
    Constants = importdata(path+"/Parameters.mat");

    Parameters.tau = Constants.tau;     %  s; Synaptic Time Delay Constant.
    Parameters.G_L = Constants.G_L;     % pS; Membrane Conductance.
    Parameters.g = Constants.g;         % pS; Synapse/Gap Conductance.
    Parameters.V_L = Constants.V_L;     % mV; Leakage Potential.
%     Parameters.L_gap = Network_Parameters.Worm279_ejunct_laplacian;
%     Parameters.W_syn = Network_Parameters.Worm279_synapse_matrix_dir;
    Parameters.L_gap = Network_Parameters.Laplacian;
    Parameters.W_syn = Network_Parameters.Synapse;
    Parameters.E = Reversal_Potential.E;
    Parameters.a_r = Constants.a_r;     %  s; Time Rise Constant.
    Parameters.a_d = Constants.a_d;     %  s; Time Delay Constant.
    Parameters.beta = Constants.beta;   % mV; Sigmoid Width.

    %% Find Threshold Voltage
    %Perturbed_Neurons = [277,279]; % The Neurons user wants to perturb.

    %I_amplitude = [2e4,2e4]; % The amplitude of the applied current. Normalized by g.

    I = zeros(N,1); % The Applied Current.
    %I(Perturbed_Neurons) = I_amplitude; 
    %I = ones(N,1).*55.798;
    
    
    I(Perturbed_Neurons) = I(Perturbed_Neurons) + I_amplitude;
    Equilibrium = Threshold_Voltage(I,path);
    

    S_eq = Equilibrium.S_eq;   % Equilibrium Synaptic Activity Variable.
    V_eq = Equilibrium.V_eq;
    V_th = Equilibrium.V_th;
    Parameters.V_th = V_th;   % mV; Threshold Voltage.

    Results.S_eq = S_eq;
    Results.V_eq = V_eq;
    Results.V_th = V_th;
    
    %% Apply Noise
    noise_var = 0.00005;
    noise_V=noise_var * randn(2*N,1);
    %% Euler Integration Set-up
    %%% Times
    start_time = 0; % 0 seconds.
    %end_time = 10; % 10 seconds.
    dt = 1e-4; % Time step-size.
    nt = round((end_time - start_time)/dt); % Total number of integration steps.
    T = start_time:dt:(end_time-dt);
    T = T'; % Time vector.

    %%% Initial Conditions
    y = zeros(2*N,1);
%     y(1:N) = V_eq;
%     y(N+1:2*N) = S_eq;
%     y(Perturbed_Neurons(1)+N) = S_eq(Perturbed_Neurons(1))+perturb; % Peturbation from equilibrium.
    y(1:N) = V_eq + perturb;
    y(N+1:2*N) = S_eq + perturb;

    %%% Storage Variable
    z = zeros(2*N,nt/(5*40));
    times = zeros(nt/(5*40),1);
    %% Euler Integrate
    j = 1;
    for i = 1:nt
        t = dt*i;
        dy = Conductance_ODE(y,N,I,Parameters);
        
        y = y + dt*dy;
        
        
        if i > int32(4/5*nt) && mod(i,40) == 0
        %if mod(i,50) == 0
            z(:,j) = y;
            times(j) = t;
            j = j + 1;
        end
    end
    V = z(1:N,:);
    S = z(N+1:2*N,:);
    Results.V = V;
    Results.S = S;
    Results.time = times;
   
end

