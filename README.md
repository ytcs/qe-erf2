# qe-erf2

## Instructions
1. Perform scf calculation using QE
2. Perform bands calculation with the same kpts using QE > bands.out
3. Clean-up bands.out by removing messages before & after the actual energy values
4. You might need to add space in front of negative values due to QE formatting
5. Run pw_exports.x (available on QE 6.2 or older) 
6. Run pp_bands.py followed by pp_ui.py
7. Run datafmt.py
8. Compile dm_{f2/f2_anisotropic/jdos}.cpp > f2.x 
9. Run f2.x {ik}
