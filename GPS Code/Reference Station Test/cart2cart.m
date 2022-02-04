function Q = cart2cart(P,t,omega,m)
% Transformation between cartesian coordinate systems from P to Q
%
% Given:
%     P       Cartesian position vector [x, y, z] in meters
%     t       Translation [dX, dY, dZ] in meters
%     omega   Rotation [omegax, omegay, omegaz] in radians
%     m       Scale [m] unitless
%     
% Returned: 
%     Q       Cartesian position vector [x, y, z]
%     
% References:
%   https://eng.mil.ru/files/PZ-90.11_final-v8.pdf

Rxyz =  [1  omega(3) -omega(2); ...
         -omega(3) 1  omega(1); ...
         omega(2) -omega(1) 1]; ...
         
Q = (1 + m)*Rxyz*P + t;