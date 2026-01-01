# Quantum-State-Transfer
Quantum State Transfer: Made with the UC Davis COSMOS program

## How to use:

Compile with `gcc -o qst.o qst.c linalg.c -lm -lblas -llapack` followed by `./qst.o`
Printed values will show the initial set up values. These can be changed via the linalg.c function "fmatrix_like_richard" referring to Professor Richard Scalettar

View results in qst-out.txt, although output.py uses matplotlib for a more "visual" method of viewing results

## What does this do?
It solves the time dependent Schrodinger in one dimension, recording probabilities of finding the given particle (via a Hermitian matrix) at any given point in the line.

## How does it do this?
We set up a hermitian matrix and provide it to this code in order to start it. Then, the code performs linear algebra operations each time step and puts the results for each step into the output file.
