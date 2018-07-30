# General_Dumbbell
Simulates a general dumbbell using several options and integrators

Core code is in Integrate.f90, including the main integrator
Generate\_Initial\_Conditions.f90 does what it says on the box
Dumbbell\_Validation\_Tests.f90 contains several stochastic 'unit' tests
Dumbbell\_util.f90 contains the random number generators, some general utility subroutines and subroutines for measuring properties of dumbbells in shear flow (given a pre-computed dumbbell ensemble)

main.f90 requires inputparameters.inp, options.inp and timestepdata.inp

inputparameters.inp has the form (sr alpha h* Q0)
timestepdata.inp has the form (No.Traj No.dtWidths output_delay; dtW1; dtW2;... ; dtWN)
options.inp has the form (VarianceReductionOpt DistributionOpt)

sr, alpha, h* and Q0 are described in my paper

No.traj is the number of dumbbells which should be in the ensemble
No.dtWidths is the number of timestep widths to be used
output_delay is how many relaxation times between outputs
dtW1,2,3,....,N must be on seperate lines and be the No.dtWidths timestep widths
VarianceReductionOpt should be set to 0 for no VR and 1 for VR
DistributionOpt prints the entire ensemble on the final timestep if it's set to 1

