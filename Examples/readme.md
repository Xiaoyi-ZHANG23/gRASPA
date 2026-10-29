# Quick Start
After compiling the code and having the executable in ```src_clean/``` folder, run the command below:
```
python run_designated_folders.py
```
Wait and see the result!
# Pytest
:memo: If you have pytest installed, after running the previous command, do:
```
pytest -s
```
You can see the energy drift and time elapsed for the test simulations you just run!

# Identity-swap regression check
```IdentitySwap-IdealGas/``` is a two-component ideal gas (all Lennard-Jones epsilons zero, no charges, empty box) whose
grand-canonical loadings are known exactly: ```<N_c> = beta * f_c * V```. It exercises the identity-swap move
(```IdentityChangeProbability```) together with insertions and deletions. Run it and compare the ```Overall: Average```
loadings with the values in ```IdentitySwap-IdealGas/README.md```, or let ```pytest -s``` do it
(```test_identity_swap_ideal_gas.py``` skips when the example has not been run).
