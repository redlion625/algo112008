function [Y,M,D]=jd2greg(jd)
% jd2greg - doesn't work
L=jd+68569;
N=4*L/146097;
L=L-(146097*N+3)/4;
I=(4000*(L+1))/1461001;
L=L-(1461*I)/4+31;
J=(80*L)/2447;
D=L-(2447*J)/80;
L=J/11;
M=J+2-12*L;
Y=100*(N-49)+I+L;

