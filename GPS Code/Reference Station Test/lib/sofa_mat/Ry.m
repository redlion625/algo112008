function r = Ry(theta,R)
s = sin(theta);
c = cos(theta);

a00 = c*R(1,1) - s*R(3,1);
a01 = c*R(1,2) - s*R(3,2);
a02 = c*R(1,3) - s*R(3,3);
a20 = s*R(1,1) + c*R(3,1);
a21 = s*R(1,2) + c*R(3,2);
a22 = s*R(1,3) + c*R(3,3);

R(1,1) = a00;
R(1,2) = a01;
R(1,3) = a02;
R(3,1) = a20;
R(3,2) = a21;
R(3,3) = a22;
r=R;