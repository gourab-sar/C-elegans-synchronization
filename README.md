# Synchronization properties in C. elegans: Relating behavioral circuits to structural and functional neuronal connectivity

## C_elegans_data

The files are arranged alphabetically
- `C_Elegans_Network_Data.mat` - Weighted neuronal connectivity data
- `Neuron_Types.mat` - Classifies the neurons by their types (sensory, interneuron, motor, polymodal)
- `Parameters.mat` - Parameters and their values used in the single compartment membrane model
- `README Worm279dir database.docx` - Description of the C. elegans data
- `Reversal_Potential.mat` - Reversal potential of the neurons
- `Worm279 Description.xlsx` - Defines neurons as excitatory/inhibitory
- `Worm279dir.mat` - Contains the binary connectivity data along with other information like birthtime, positions, labels

We refer the readers to the papers [Structural Properties of the Caenorhabditis elegans Neuronal Network](https://doi.org/10.1371/journal.pcbi.1001066) for description of the C. elegans structural connectome, and to [Low-dimensional functionality of complex network dynamics:Neurosensory integration in the Caenorhabditis elegans connectome](https://doi.org/10.1103/PhysRevE.89.052805) for the description of the single compartment membrane model and associated parameters.

## Main_matlab_codes
Will update here...

## Python_codes
Will update here...

## Supporting_matlab_functions
We recommend not making changes to these auxiliary functions

### Supporting functions for the main model
- `Conductance_ODE.m` - Defines the conductance based single compartment membrane model equations
- `Find_Hopf_Bifurcation_Stepper.m` - Sets up the Hopf bifurcation method
- `Run_Model.m` - Sets up the Euler-Maruyama integration method. Change step size, noise amplitude, initial conditions, specifications about data to be stored here.
- `Threshold_Voltage.m` - Calculates the threshold voltages of the neurons. 

### Supporting functions for the Weighted Stochastic Block Model (WSBM) method
- `Adj2Edg.m`, `biwsbm.m`, `calc_logEvidence.m`, `Edg2Adj.m`, `setup_distr.m`, `wsbm.m`
  
We refer the readers to [The Weighted Stochastic Block Model](https://aaronclauset.github.io/wsbm/) by Aaron Clauset for the documentation and proper implementation of WSBM

### Supporting functions for the generalized Louvain method
- `genlouvain.m`, `iterated_genlouvain.m`
  
We refer the readers to [A generalized Louvain method for community detection implemented in MATLAB](https://github.com/GenLouvain/GenLouvain) for the documentation and proper implementation of the generalized Louvain method

**Important:** We strongly recommend to put all the data from the *C_elegans_data* folder, and all the codes and supporting functions from the *Main_matlab_codes*. *Python_codes*, and *Supporting_matlab_functions* folders into a single folder before running any simulations.

