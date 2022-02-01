function ts = emissionTime(R,tr,te,a0,a1,a2)
% emissionTime
% Given:
%     R           pseudorange
%     tr          time of receiver
%     a0,a1,a2    2nd order SV clock offset polynomial coefficients
% 
% Returned:
%     ts          Time of signal emission from satellite
%
% Reference:
%   https://gssc.esa.int/navipedia/index.php/Emission_Time_Computation

ts = tr - R/Const.CMPS;
ts = ts - svClockOffset(ts,te,a0,a1,a2);
    