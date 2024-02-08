%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function [minc,volfin,gray] = SurrOpt_top88_inner_Fin(nelx,nely,volfrac,penal,rmin,ft)

%% INITIALIZE 
addpath('MMA')
nele = nelx*nely;
%% MATERIAL PROPERTIES
E0 = 2;
Emin = 1e-9;
nu = 0.3;

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
% Modification 1 - Defining betaH and filter densities before start of loop
betaH = 5;
eta = 0.5;
if ft == 1 || ft == 2
    xPhys = x;
elseif ft == 3
    xTilde = x;
    xPhys = (tanh(betaH*eta)+tanh(betaH*(xTilde-eta)))./...
    (tanh(betaH*eta)+tanh(betaH*(1-eta))); %250 line Heaviside
end
loop = 0;
change = 1;

%% INITIALIZE MMA OPTIMIZER
m     = 1;                % The number of general constraints.
n     = nele;             % The number of design variables x_j.
xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);             % xval, one iteration ago (provided that iter>1).
xold2 = x(:);             % xval, two iterations ago (provided that iter>2).
low   = ones(n,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

%% START ITERATION
while change > 0.001
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  % Modification 2
  elseif ft == 3
    dx = betaH*(1-tanh(betaH*(xTilde-eta)).^2)./...
    (tanh(betaH*eta)+tanh(betaH*(1-eta))); %250 Line Heaviside
    dc(:) = H*(dc(:).*dx(:)./Hs);
    dv(:) = H*(dv(:).*dx(:)./Hs);
  end
  
    %% METHOD OF MOVING ASYMPTOTES
    xval  = x(:);
    f0val = c;
    df0dx = dc(:);
    fval  = sum(xPhys(:))/(volfrac*nele) - 1;
    dfdx  = dv(:)' / (volfrac*nele);
    [xmma, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
    mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d,betaH);
    % Update MMA Variables
    xnew     = reshape(xmma,nely,nelx);
    xTilde(:) = (H*xnew(:))./Hs;
    xPhys = (tanh(betaH*eta)+tanh(betaH*(xTilde-eta)))./...
    (tanh(betaH*eta)+tanh(betaH*(1-eta))); %250 line Heaviside
   
    xold2    = xold1(:);
    xold1    = x(:);
  
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
%   tiledlayout(3,1)
%   nexttile
%   %Design variable
%   colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;
%   nexttile
%   %Filter
%   colormap(gray); imagesc(1-xTilde); caxis([0 1]); axis equal; axis off; drawnow;
%   nexttile
  %Heaviside/Physical
  figure(2)
  colormap(jet); imagesc(xPhys); caxis([0 1]); axis equal; axis tight; axis off; drawnow;
%   pause(1);
  %Modification 4 - continuation scheme for the regularization parameter b
%   if ft == 3 && betaH<512 &&...
%           (loopbeta >= 50 || change <= 0.01)
%       betaH = 2*betaH;
%       loopbeta = 0;
%       change = 1;
%       fprintf('Parameter betaH increased to %g.\n',betaH);
%   end 

    %===========APPLYING CONTINUATION=============%
    if loop >= 600 || change <= 0.001
        break
    end

end

minc = c;
volfin = mean(xPhys(:));
gray = sum(sum(xPhys<=0.95 & xPhys>=0.05)) / nele; % percent of grey / total ele

figure(2)
end 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

