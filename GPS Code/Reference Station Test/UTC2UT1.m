function UT1 = UTC2UT1(UTC,UT1UTC)
UT1=UTC+UT1UTC/86400;
% Disclaimer: This software is based off of SOFA c source code using the
% same definitions but is a derivative work and differs from the SOFA
% source code with the known differences:
% 1) IAU resolutions implemented in MATLAB instead of C
% 2) Less accurate methods of conversion between time scales (TT,UT1,UTC)