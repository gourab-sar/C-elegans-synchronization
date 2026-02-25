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

## Supporting_matlab_functions
We recommend not making changes to these auxiliary functions

### Supporting functions for the main model
- `Conductance_ODE.m` - Defines the conductance based single compartment membrane model equations
- `Find_Hopf_Bifurcation_Stepper.m` - Sets up the Hopf bifurcation method
- `Run_Model.m` - Sets up the Euler-Maruyama integration method. Change step size, noise amplitude, initial conditions, specifications about data to be stored here.
- `Threshold_Voltage.m` - Calculates the threshold voltages of the neurons. 
