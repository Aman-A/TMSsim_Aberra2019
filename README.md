# README for TMSsim_Aberra2019

Simulation code for Aberra, A., Wang, B., Grill, W. M., and Peterchev, A. V. (2019). Simulation of transcranial magnetic stimulation in head model with morphologically-realistic cortical neurons. Brain Stimul. In press.

Dependencies:
- NEURON v.7.4 or later (may run on some older versions)
- MATLAB r2016 or later (should run on most older versions)
    Certain functions for generating figures require additional MATLAB toolboxes
        - Curve Fitting Toolbox: required by `makeSlices.m` (used in `plotFig2b.m`, `plotFig4.m`, `plotFig5.m`)    

To run NEURON simulations of single neurons:

1) Start by compiling NEURON mechanisms in (./nrn/mechanisms/) folder, then drag `special`/`nrnmech.dll` file (Unix/Windows) into (./nrn/) folder

On MacOS/Linux:
Execute following code to make executable:
`chmod +x mosinit.command`

2) Initialize cell data by running `InitializeCellData.m` in the (./mat/) folder
    Upon completion, this should generate cell_data/cell_models.mat, as well as cell_data/maxH and cell_data/maxHlin folders with following sub-folders: areas/, branchorders/, comp_types/, coordinates/, diams/, parent_inds/, secnames/, sectypes/. Each folder contains the corresponding data for the 25 cell models included in the E-field coupled population model
    Also generates 6 `NeuronPop` structs for the "maxH" model version (`nrn_model_ver`), stored in "nrn_pop[pop_ind]_maxH.mat", where pop_ind is 1–6, and 1 `NeuronPop` for the "maxHlin" model version (L2/3 pyramidal cell with straight axon)
3) Interpolate E-field FEM solutions on compartments of neuron model populations to simulate thresholds by running `interpEfieldsAll.m` in the (./mat/) folder.
    This should generate arrays of E-field vectors for each neuron in the neuron populations specified in `pop_inds` (1–6). 
    WARNING: this step may take ~2 hours if run on a single CPU serially. If run on multi-core machine, the function interpEfield will automatically parallelize the computation across all available CPUs.
    Additionally, E-field data takes up 50.3 MB per neuron, summing 1.25 GB for all 25 neurons in each population. Interpolating E-field vectors for all 6 populations, 2 E-field directions (for Medtronic MC B70 coil), and 2 E-field solutions (M1_PA_MCB70 and M1_PA_Magstim70mm) requires ~22 GB, so it is not recommended to generate all E-field, unless simulations with all cell types/populations/E-fields are going to be run

4) Once cell_data and nrn_efields data has been generated, individual neural thresholds with TMS-induced E-field and waveforms can be run using set_test_params.m script to save the TMS waveform (tvec.txt/Evec.txt) and E-field (Er.txt) distribution to (./nrn/params/test) to be used in NEURON simulation, and specifying parameters in (./params/test/defineParams.hoc). Default settings initialize a L23 pyramidal cell for uniform E-field stimulation (`load_potentials = 0`). Set `load_potentials` = 1 to use saved Er.txt. NOTE: Er.txt must be saved for correct neuron model, which is designated by `cell_id` (L23 pyramidal cells are 6-10).

The NEURON simulation can then be run by running (./nrn/.mosinit.command) from the terminal (MacOS/Linux) or 


# SLURM
Example code for parallelizing simulations on linux cluster with SLURM scheduler is in slurm/ folder

