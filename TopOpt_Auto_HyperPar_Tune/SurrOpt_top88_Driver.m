%============THIS PROGRAM DEMO THE STRUCTURE TO AUTOMATICALLY TUNE THE
%HYPERPARAMETERS OF TOPOLOGY OPTIMIZATION ALGORITHMS=============%

close all
clear
clc

%INITIALIZATION - Input constant hyperparameters for the inner topology optimization (TO) algorithm 
%The TO algorithm here is based on the 88-line MATLAB code (Andreassen et
%al., 2011)
nelx = 150; %length of the domain (in elements)
nely = 50; %height of the domain (in elements) 
volfrac = 0.5; %prescribed volume fraction
rmin = 10/3; %Filtering radius 
ft = 3; %Filtering scheme - 3 = Heaviside 

%Provide the bounds for the hyperparamaters to tune 
LB = [1];
UB = [10];

%Provide the settings for Surrogate Optimization before execution
minpoints = 2;
maxfunceval = 10;

%Provide other necessary settings - refer to MATLAB documentation
options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
    'ConstraintTolerance',1e-3,'MaxFunctionEvaluations',maxfunceval,'MinSurrogatePoints',minpoints,...
    'MinSampleDistance',0.1,'Display','iter');

%Call the surrogateopt function in MATLAB
[A,f,exitflag,output,trials] = surrogateopt(@(A)obj(A,nelx,nely,volfrac,rmin,ft),LB,UB,[],[],[],[],[],options);
%After Surrogate Optimization has finished, the inner TO algorithm is
%called one more time to plot and acquire the outer objective components
[minc,volfin,gray] = SurrOpt_top88_inner_Fin(nelx,nely,volfrac,A(1),rmin,ft);
outerobj(:,1) = [minc;volfin;gray];

%Save outputs and figures
save('A.mat',"A"); %A is the vector of tuned hyperparameters
save('outerobj.mat',"outerobj"); %outerobj is the components of the outer objective function
saveas(figure(2),'TopFin.fig'); %Save the plotted topology as .fig file

function func = obj(A,nelx,nely,volfrac,rmin,ft)
%This function is called at each surrogate point, which performs a full TO
%and outputs the components for the objective function, namely: (1) minc =
%compliance; (2) volfin = volume fraction of the design; (3) gray =
%grayness metric of the design 
[minc,volfin,gray] = SurrOpt_top88_inner(nelx,nely,volfrac,A(1),rmin,ft);

%The variable that determines the weighting of component (2) and (3), written as zeta in the paper
fac = 100; 

%The objective components are combined to form the outer objective
func.Fval  = minc + fac*(volfin - volfrac)^2 + fac*gray; 
end