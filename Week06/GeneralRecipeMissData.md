# General recipe for missing data problem

For future reference, I provide the general recipe for parameter estimation in missing data problems using the EM algorithm. Here is the recipe:

1. Determine the full data log-likelihood. This function should be analytically tractable.  
2. Determine the Q-function, i.e. calculate the E-step.   
3. Maximize the Q function, i.e. carry out the M-step.
4. Implement the iterative procedure.
5. Decide on initial parameter values and a stopping criteria.
6. Check that the data log-likelihood increases in every iteration.  
7. Running the algorithm multiple times with different initial parameters is also useful because the algorithm only converges to a local maximum.
