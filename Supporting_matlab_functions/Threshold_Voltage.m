%% Find Threshold Voltage.
function Results = Threshold_Voltage(I,path)
    %%% Finds the threshold voltage used for the sigmoid function phi.
    %   The vector I is the applied current.
    
    N = 279; % Total Number of Neurons.
    Parameters = importdata(path+"/Parameters.mat");
%     Network_Data = importdata(path+"/Worm279dir.mat");
    Network_Data = importdata("C_Elegans_Network_Data.mat");
    Reversal_Potential = importdata(path+"/Reversal_Potential.mat");

    %%% Voltage ODE Parameters.
    G_L = Parameters.G_L;
    g = Parameters.g;
    tau = Parameters.tau;
    V_L = Parameters.V_L;
%     L = Network_Data.Worm279_ejunct_laplacian;
%     W_syn = Network_Data.Worm279_synapse_matrix_dir;
    L = Network_Data.Laplacian;
    W_syn = Network_Data.Synapse;
    E = Reversal_Potential.E;
    
    %%% Synapse ODE Parameters.
    a_r = Parameters.a_r;
    a_d = Parameters.a_d;
    beta = Parameters.beta;
    phi_eq = 0.5;

    %%% Solve for the Equilibrium Synaptic Activity Variable, S_eq.
    S_eq = ones(N,1).*(phi_eq/((a_d/a_r) + phi_eq));
   
    %I = zeros(N,1);
    %I([146,147]) = 2e4;

    %%% Solve for the Equilibrium Potential, V_eq.
    b = I + (G_L/g)*V_L + W_syn'*(S_eq.*E);
    A = (G_L/g)*eye(N) + diag(W_syn'*S_eq)+ L;
%     b = I + ((G_L/g)*V_L + W_syn'*(S_eq.*E))/tau;
%     A = ((G_L/g)*eye(N) + diag(W_syn'*S_eq)+ L)/tau;
    V_eq = A^-1*b;

    %%% Solve for the Threshold Potential, V_th.
    V_th = V_eq; % + 1/beta*log(1/phi_eq - 1);
    
    Results.S_eq = S_eq;
    Results.V_eq = V_eq;
    Results.V_th = V_th;
end