# Human-Informed Topology Optimization 
This is the public code that was featured in the paper "Human Infored Topology Optimization: Interactive Application of Feature Size Controls" (Ha, D. and Carstensen, J.V., 2022). 
The program seeks to perform topology optimization of a 3-point bending MBB beam, where the user is able to interact with the design during the optimization process.
The following folder and file are provided and required to run the program:
1. Folder ``` MMA```
2. File ``` hiTop_MBB_1ROI_Fin.m```
Once the files are downloaded, the program can be run on MATLAB with the command to call the function:
```
hiTop_MBB_1ROI_Fin(nelx,nely,volfrac,penal,rmin,betaH,maxit)
```
Where ```nelx, nely``` define the size of the design domain, ```volfrac``` defines the maximum amount of materials that are allowed, as a fraction of the total space, ```penal``` defines the penalty factor in SIMP, ```betaH``` defines the steepness factor beta in the Heaviside Projection Method, and  ```maxit``` defines the maximum number of iteration that the code is allowed to run.
As a starting point, the user is recommended to try the following:
```
hiTop_MBB_1ROI_Fin(240,80,0.5,3,3.2,25,700)
```
Then, the program will run automatically for the first 50 iterations, where the features begin to form. 
![50th iteration of the design](https://i.imgur.com/woUrmqN.png)

At the 50th iteration, the program will stop, and focuses on the Figure window of the design. Then, the user will be asked to highlight one (1) region of interest (ROI) where the user believes a member could be too thin. The cursor on the Figure window will switch to a crosshair shape. The user can then draw an elliptical shape that surrounds a feature of concern, as shown:
![HiTop ROI specification](https://user-images.githubusercontent.com/112650617/218153118-6cc9ddf5-de6c-484e-b90f-184c038d3140.png)
Once the user finishes drawing the ROI, the shape is finalized by double clicking the left mouse button on the center of the ROI. Once this is done, in the Command Window of MATLAB, the user is prompted to enter the following: 
```
Enter the number of contour lines
```
This input indicates how many contours, generating inwards, will the feature size map change over. The user type in a number, and press Enter. For the demonstration purpose, the user is recommended to enter 3.
Next, the user is prompted to enter the following:
```
Input rmin of ROI
```
This input indicates the minimum feature size within the ROI. Please note that this input has to be bigger than the overall feature size of the design that was specified right at the beginning with the variable ```rmin```. In this case, the user is recommended to enter 6.4.
Finally, the user is prompted to enter the following:
```
Input lambda between domain and ROI
```
This input defines the degree of transition that the feature size outside can transition into the ROI. A higher number indicates a gradual transition. Conversely, a lower number indicates a sharper transition, and the user will notice more dramatic changes in the response. For the example case, the recommended value is 3.
Then, the program resumes, and the new topology with the human input becomes:
![image](https://user-images.githubusercontent.com/112650617/218154272-e702fb65-e398-45f1-8c73-f66170545e88.png)
