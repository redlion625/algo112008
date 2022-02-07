function [X,Y,Z] = refSysTransform(X0,Y0,Z0,dX,dY,dZ,omegaX,omegaY,omegaZ,m)
% dX,dY,dZ - translation/offset [m]
% omegaX,omegaY,omegaZ - rotation [rad]
% m - scale (fractional)

% Notes: reverse signs of all transformation elements for inverse transformation
P = [X0;Y0;Z0];
Q=(1+m)*[1 omegaZ -omegaY;-omegaZ 1 omegaX;omegaY -omegaX 1]*P+[dX;dY;dZ];
X=Q(1);Y=Q(2);Z=Q(3);


