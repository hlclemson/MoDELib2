# note
- This script determines the critical resolved shear stress (CRSS) of the AlMg5 alloy system using a bisection-search algorithm under correlated analytical solid solution noise and correlated stacking fault noise. First, the system is minimized for the number of steps specified by the output frequency. An initial external stress of 1 MPa is then applied; if the dislocations do not move, the stress is increased in 100 MPa increments.
- To set how many random seeds to test, edit the variable noise_seed_to_test in main.py.
