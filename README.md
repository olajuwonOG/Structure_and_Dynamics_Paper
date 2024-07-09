# Structure_and_Dynamics_Paper
Lammps and python script for structure and dynamics of dispersed polymer melts


- mkinput_monodisperse.py: python code to generate initial configurations of monodisperse melt
- mkinput_monodisperse.py: python script to generate initial configurations of polydispersed melt
- in.poly: lammps input file for equilibration and production
- msd.py: python code for analyzing dynamics (monomer and center of mass mean squared displacement for all chains and test chain)
- e2etest.py: python script to call -pppmd.py for computing end-to-end vector autocorrelation (E2E ACF)
- pppmd.py: post processing script for calculation - used for E2E ACF in this case.
- radiusOfGyrationPolydispersity.py: python script for calculating average radius of gyration of all chains and test chains
- Dynamic_Structure_Factor.py: calculates dynamic structure factors of test chains in melt
- Static_Struct_Factor.py: calculates the static structure factors for test chains and bulk
