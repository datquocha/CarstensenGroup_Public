README
# Automatic hyperparameter tuning of topology optimization (TO) algorithms

This is the Demo of the approach presented in the study "Automatic hyperparameter tuning of topology optimization algorithms using surrogate optimization". 

Here, the hyperparameters of the 88-line MATLAB code (Andreassen et al., 2011), solving the MBB compliance minimization problem, is automatically tuned. 
In this Demo, the SIMP penalty $p$ is optimized to acquire a topology solution that (1) minimizes compliance; (2) has crisp feature boundaries; (3) adhere to the prescribed volume fraction constraint. 

The Surrogate Optimization module in MATLAB is used to solve the nested optimization problem. The inner problem is the standard TO problem, framed by the outer problem that has the hyperparameter as the design variable that satisfy requirements (1), (2), (3) above. 

The following folder and file are provided and required to run the program:

1. Folder ```MMA```
2. File ```SurrOpt_top88_Driver.m```
3. File ```SurrOpt_top88_inner.m```
4. File ```SurrOpt_top88_inner_Fin.m```

Folder 1 and files 2-4 should be downloaded and stored in one single directory to ensure proper execution.

The file ```SurrOpt_top88_Driver.m``` is the main program, which in turn calls ```SurrOpt_top88_inner.m``` and ```SurrOpt_top88_inner_Fin.m```. 
The program can be readily run to tune the hyperparameter SIMP penalty $p$. 
The code is structured to present the applicability of the framework to a wide range of hyperparameters, TO programs, and outer objective function. 

We welcome any questions and comments. Please send them to the email ```datha@mit.edu```.

Disclaimer: The authors reserves all rights but do not guaranty that the code is free from errors. Furthermore, we shall not be liable in any event caused by the use of the program. 
