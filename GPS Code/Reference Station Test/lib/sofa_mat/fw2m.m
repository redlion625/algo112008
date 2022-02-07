function r = fw2m(gamb,phib,psi,eps)

r=eye(3);
r=Rz(gamb,r);
r=Rx(phib,r);
r=Rz(-psi,r);
r=Rx(-eps,r);