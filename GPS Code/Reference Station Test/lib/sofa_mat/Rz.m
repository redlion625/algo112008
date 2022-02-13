function r = Rz(psi,R)

s = sin(psi);
c = cos(psi);

a00 =   c*R(1,1) + s*R(2,1);
a01 =   c*R(1,2) + s*R(2,2);
a02 =   c*R(1,3) + s*R(2,3);
a10 = - s*R(1,1) + c*R(2,1);
a11 = - s*R(1,2) + c*R(2,2);
a12 = - s*R(1,3) + c*R(2,3);

R(1,1) = a00;
R(1,2) = a01;
R(1,3) = a02;
R(2,1) = a10;
R(2,2) = a11;
R(2,3) = a12;
r=R;