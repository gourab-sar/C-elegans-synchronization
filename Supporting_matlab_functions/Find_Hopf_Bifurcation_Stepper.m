%% Finding the Hopf Bifurcation point for the applied current
function EV = Find_Hopf_Bifurcation_Stepper(I_test,Perturbed_Neurons,path)
    %Perturbed_Neurons = [277,279];

    %%% Loading Relevant Data.
    N = 279; % Total Number of Neurons.
%     Network_Parameters = importdata(path+"/Worm279dir.mat");
    Network_Parameters = importdata("C_Elegans_Network_Data.mat");
    Reversal_Potential = importdata(path+"/Reversal_Potential.mat");
    Constants = importdata(path+"/Parameters.mat");

    %%% Parameters.
    Parameters.tau = Constants.tau;     %  s; Synaptic Time Delay Constant.
    Parameters.G_L = Constants.G_L;     % pS; Membrane Conductance.
    Parameters.g = Constants.g;         % pS; Synapse/Gap Conductance.
    Parameters.V_L = Constants.V_L;     % mV; Leakage Potential.
%     Parameters.L_gap = Network_Parameters.Worm279_ejunct_laplacian;         % Laplacian Matrix.
%     Parameters.W_syn = Network_Parameters.Worm279_synapse_matrix_dir; % Synaptic Adjacency Matrix.
    Parameters.L_gap = Network_Parameters.Laplacian;         % Laplacian Matrix.
    Parameters.W_syn = Network_Parameters.Synapse; % Synaptic Adjacency Matrix.
    Parameters.E = Reversal_Potential.E;         % Reversal Potentials.
    Parameters.a_r = Constants.a_r;     %  s; Time Rise Constant.
    Parameters.a_d = Constants.a_d;     %  s; Time Delay Constant.
    Parameters.beta = Constants.beta;   % mV; Sigmoid Width
    
    %%% Find Threshold Voltage.
    I = zeros(N,1); % The Applied Current.
    I(Perturbed_Neurons) = I_test;

    %Equilibrium = Run_Model(Perturbed_Neurons,I_amplitude,Plotting);
    Equilibrium = Threshold_Voltage(I,path);
    
    
    S_eq = Equilibrium.S_eq;   % Equilibrium Synaptic Activity Variable.
    V_eq = Equilibrium.V_eq;
    V_th = Equilibrium.V_th;
    Parameters.V_th = V_th;

    %%% Initial Conditions
    x_0 = zeros(2*N,1);
    x_0(1:N) = V_eq;
    x_0(N+1:2*N) = S_eq;
%     x_0(Perturbed_Neurons+N,1) = S_eq(Perturbed_Neurons,1)+2.5e-3; % Peturbation from equilibrium.

    %%% Find EigenVectors.
    options = optimset('Display','off');

    fun = @(x)Conductance_ODE(x,N,I,Parameters);

    [~,~,exitflag,~,Jacob] = fsolve(fun,x_0,options);
    Eigen_Values = eig(Jacob);

    if exitflag < 0
        return
    end

    EV = Eigen_Values;
    
end