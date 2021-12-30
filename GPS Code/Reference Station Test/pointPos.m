function [XYZ,v,DOP] = pointPos(L_PR,Xs,Ys,Zs,dTsv,xapprox,yapprox,zapprox)
% pointPos - calculate least squares solution for GNSS point positioning
% Inputs:
%   L_PR - n x 1 vector of pseudoranges
%   Xs,Ys,Zs - n x 1 satellite positions
%   dTsv - n x 1 satellite clock offset
%   xapprox,yapprox,zapprox - initial estimate of receiver position
% Outputs:
%   XYZ - 3 x 1 vector of positions in WGS84 cartesian coordinates
%   v - n x 1 residuals for each pseudorange
%   DOP - 4 x 1 dilution of precision (LS error) of x,y,z and dTr

% Notes: DOP not scaled by a posteriori  factor

c = 2.99792458e8;

%The initial approximate of the user satellite position for the observation file
%xapprox = 0;  % Obs_data.pos_xyz(1,1);
%yapprox = 0;   %Obs_data.pos_xyz(1,2);
%zapprox = 0;   %Obs_data.pos_xyz(1,3);

% Xs = Observations{i,2}(:,14); %The current x Earth-fixed coordinates of SV antenna phase center
% Ys = Observations{i,2}(:,15); %The current y Earth-fixed coordinates of SV antenna phase center
% Zs = Observations{i,2}(:,16); %The current z Earth-fixed coordinates of SV antenna phase center

L=L_PR;
%Setting the initial approximations of the receiver positions and the
%receiver clock offset
xcurrent = xapprox;  % Obs_data.pos_xyz(1,1);
ycurrent = yapprox;   %Obs_data.pos_xyz(1,2);
zcurrent = zapprox;   %Obs_data.pos_xyz(1,3);
dtrcurrent = 0;
%xcurrent = xapprox;   ycurrent = yapprox;   zcurrent = zapprox;   dtrcurrent = 0;
count = 0;

n=numel(L);
PR=zeros(n,1);
%Performing a nonlinear least squares adjustment to determine the least
%squares estimates of the receiver x, y, z position and the receiver
%clock offset
while(count < 40)
    clear PR A mis
    x0 = xcurrent;   y0 = ycurrent;   z0 = zcurrent;  dtr0 = dtrcurrent; %Setting the current estimates of the unknown parameters to the current least squares estimates of the parameters
    
    for j = 1:n %A for-loop to create the first design matrix
        %Determining the geometric range
        geometricexpression = (Xs(j) - x0)^2 + (Ys(j)- y0)^2 + (Zs(j) - z0)^2;
        geometricinitial = sqrt(geometricexpression); %+c(dts-dtr)
        
        %The deriviatives for the X, Y, and Z receiver components and the
        %receiver clock offset
        derivx = -((Xs(j) - x0)/geometricinitial);
        derivy = -((Ys(j) - y0)/geometricinitial);
        derivz = -((Zs(j) - z0)/geometricinitial);
        derivdtr = c;
        
        %Populating the first design matrix
        A(j,1) = derivx;
        A(j,2) = derivy;
        A(j,3) = derivz;
        A(j,4) = derivdtr;
        
        %Determining the current uncorrected pseudorange for the
        %approximatations of the receiver position and clock offset
        PR(j,1) = sqrt(geometricexpression) + c * (dtr0 - dTsv(j));  %- dTsv(j)
    end
    
    mis = L - PR; %Determining the misclosure between the current estimate of the pseudorange and the corrected pseudorange
    %P = eye(obssize(1));
    
    %Determining the updating value (deltax) to the current
    %approximations of the receiver position and clock offset
    Inverseterm = (A')*A;
    Inverse = inv(Inverseterm);
    deltax = Inverse*(A')*mis;
    xo = [x0; y0; z0; dtr0]; %The current approximations of the receiver position and clock offset
    xhat = xo + deltax; %The current least squares estimates of the receiver position and clock offset
    
    %Setting the current approximations of the receiver position and clock offset to the current least squares estimates
    xcurrent = xhat(1);
    ycurrent = xhat(2);
    zcurrent = xhat(3);
    dtrcurrent = xhat(4);
    %disp(j);
    thresh = [0.01; 0.01; 0.01; 1e-4]; %Determining the threshold at which the least squares loop with break
    count = count + 1; %Keeping track of how many iterations the loop makes
    if(all(abs(deltax) < thresh)) %If the delta values are less than the prescribed threshold values, the loop will break
        fprintf("Loop is broken");
       % resultcheck(i) = count;
        Final = xhat; %Storing the least squares estimates of the X, Y, and Z components of the receiver's position and the receiver clock offset
        v = A*deltax - mis; %Calculating the residuals
        %Observations{i,2}(:,17) = v; %Storing the residuals
        aposteriori = ((v')*v)/(n - 4); %Determining the a-posteriori variance factor
        %             Final(i,2) = aposteriori*Inverse;
        break;
    end
end
XYZ=Final(1:3); %only point position returned, receiver offset discarded
% DOP = Inverseterm; %Determining the DOP values
DOP = diag(Inverse); % not sure why scaling by posteriori 