%%%% A MODIFIED 88 LINE TOPOLOGY OPTIMIZATION CODE - hiTop - February, 2023 %%%%
%%Authors: Dat Ha (datha@mit.edu), Josephine Carstensen (jvcar@mit.edu)
%%Carstensen Group, MIT CEE 

%Example Command: hiTop_MBB_1ROI(240,80,0.5,3,3.2,25,500)
%ROI - contour = 3, rmin2 = 6.4, lambda = 3
function hiTop_MBB_1ROI(nelx,nely,volfrac,penal,rmin,betaH,maxit)

close all
clc

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
% iH = ones(nelx*nely*(2*(ceil(max(rmin1,rmin2))-1)+1)^2, 1); %ceil(max(rmin1,rmin2))
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    r(j1,i1) = rmin;
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

figure(1)
colormap('summer');
tiledlayout(1,2);
nexttile;
imagesc(r); title(['Original rmin']);
axis equal; axis tight; axis off; colorbar('Location','southoutside'); 

%% INITIALIZE ITERATION
x = repmat(0.01,nely,nelx);
xTilde = x;
xPhys = 1 - exp(-betaH*xTilde) + xTilde*exp(-betaH);

loop = 0;
loopbeta = 0;
change = 1;
mark = 0;

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
  loopbeta = loopbeta + 1;

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
  dx = betaH*exp(-betaH*xTilde) + exp(-betaH);
  dc(:) = H*(dc(:).*dx(:)./Hs);
  dv(:) = H*(dv(:).*dx(:)./Hs);
  
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
    xPhys = 1 - exp(-betaH*xTilde) + xTilde*exp(-betaH);
    
    xold2    = xold1(:);
    xold1    = x(:);
  
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
  
  %% PRINT RESULTS
  figure(2)
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);

  %% PLOT DENSITIES
  colormap(jet); imagesc(xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  if loop == 50 && mark == 0
      figure(2)
      %-----------------ROI 1--------------%
      fprintf('Arrived at loop %g.\n ',loop);
      fprintf('Highlight ROI 1 - double click on ROI when finished editing \n');
      roi = drawellipse('Color','white');
      %Wait until user finishes editing ROI
      l = addlistener(roi,'ROIClicked',@clickCallback);
      uiwait;
      delete(l);
      
      contno = input('Enter number of contour lines \n');
      
      rnew = input('Input rmin of ROI \n');
      
      lambda = input('Input lambda between domain and ROI \n');
      
      %Retrieve relevant parameters in drawn ROI
      center = roi.Center;
      axes = roi.SemiAxes;
      angle = roi.RotationAngle;
      
      %Generate arrays of x and y for inROI checks
      xcoord = zeros(1,(nelx*nely));
      ycoord = zeros(1,(nelx*nely));
      count = 1; 
      for i=1:nely
          for k=1:nelx
              xcoord(count)=k;
              count = count+1;
          end 
      end 
      count = 1;
      for i=1:nely
          for k=1:nelx
              ycoord(count)=i;
              count = count+1;
          end 
      end 

      %Generate rmin function from input lambda
      xfunc = 1:1:nelx;
      rfunc = zeros(1,nelx);
      for i=1:nelx
        %Note: Only works right now for rnew > rmin
        rfunc(i) = 0.5*(rnew+rmin) + 0.5*(rnew-rmin)*tanh((i-nelx/2)/lambda);
      end 

        for i=1:length(rfunc)
            if abs(rfunc(i)-rmin) < 0.01
                start = i;
            elseif abs(rfunc(i)-rnew) > 0.01
                fin = i;
            end 
        end 
    
        %Array of r for plotting ROI
        rroi = r; 
        d = fin - start;
        contgap = ceil(d/contno);
    
          tf = inROI(roi,xcoord,ycoord);
          for i=1:length(tf)
              if tf(i) == 1
                  r(ycoord(i),xcoord(i)) = rfunc(start);
                  rroi(ycoord(i),xcoord(i)) = rfunc(start);
              end 
          end 
          
        for z=1:contno
              %New concentric contour - Generating inwards
              roinew = drawellipse('Center',center,'SemiAxes',axes-3*z,'RotationAngle',angle,'Color','white');   
              tfnew = inROI(roinew,xcoord,ycoord);
              for i=1:length(tfnew)
                  if tfnew(i)==1 
                      r(ycoord(i),xcoord(i)) = rfunc(start + contgap*z);
                  end 
              end 
        end

    keyboard

    %% PREPARE FILTER - AFTER UPDATING NEW RMIN

    figure(1)
    nexttile;
    imagesc(r); axis equal tight;
    title(['Processed ROI with given degree of transition (lambda)']);
    colormap('summer');
    colorbar('Location','southoutside'); 
    
    figure(2)

    for i1 = 1:nelx
      for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(r(j1,i1))-1),1):min(i1+(ceil(r(j1,i1))-1),nelx)
          for j2 = max(j1-(ceil(r(j1,i1))-1),1):min(j1+(ceil(r(j1,i1))-1),nely)
            e2 = (i2-1)*nely+j2;
            k = k+1;
            iH(k) = e1;
            jH(k) = e2;
            sH(k) = max(0,r(j1,i1)-sqrt((i1-i2)^2+(j1-j2)^2));
          end
        end
      end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);

    % Reseting design variables and xPhys
    x = repmat(0.01,nely,nelx);
    xTilde = x;
    xPhys = 1 - exp(-betaH*xTilde) + xTilde*exp(-betaH);
   
    loop = 0;
    loopbeta = 0;
    change = 1;
    mark = 1;
    
  end 
 
  %Exit condition
  if loop == maxit
      xInitial = xPhys;
      save('x0.mat','xInitial');
      fprintf('Reached loop limit %g.\n ',loop);
      break
  end

  %Increase penal by 1 every 50 iterations
  if mod(loop,50) == 0 && mark == 1 && loop>0 && penal < 10
        penal = penal + 1;
        fprintf('Arrived at loop %g., Current penalty %g\n ',loop,penal);
  end 
end
end 

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end 

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

