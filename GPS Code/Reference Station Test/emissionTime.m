function ts = emissionTime(R,tr,te,a0,a1,a2)
% emissionTime
% Given:
%     R           Pseudorange
%     tr          Time of receiver in seconds of some start epoch
%     te          Time of navigation reference epoch
%     a0,a1,a2    2nd order SV clock offset polynomial coefficients
% 
% Returned:
%     ts          Time of signal emission from satellite in seconds from
%                 start epoch
%
% Notes:
%       start reference epoch of ts and tr must be the same
%
%       Relativistic correction not applied, GLONASS implements this in
%       clock bias and drift terms from navigation message
%
% Reference:
%   https://gssc.esa.int/navipedia/index.php/Emission_Time_Computation

ts = tr - R/Const.CMPS;
ts = ts - svClockOffset(ts,te,a0,a1,a2);
    