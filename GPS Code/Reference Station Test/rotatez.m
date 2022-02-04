function [x1,y1,z1] = rotatez(x0,y0,z0,S)
% rotatez(x0,y0,z0,S) - rotates 3D vector by angle S in radians about Z
%                       axis in a right hand system

x1 = x0*cos(S) - y0*sin(S);
y1 = x0*cos(S) + y0*cos(S);
z1 = z0;