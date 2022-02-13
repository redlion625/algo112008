function [gamb,phib,psib,epsa] = pfw06(date1,date2)


% Interval between fundamental date J2000.0 and given date (JC). 
t = ((date1 - Const.DJ00) + date2) / Const.DJC;

% P03 bias+precession angles. 
gamb = (    -0.052928     + ...
       (    10.556378     + ...
       (     0.4932044    + ...
       (    -0.00031238   + ...
       (    -0.000002788  + ...
       (     0.0000000260 ) ...
       * t) * t) * t) * t) * t) * Const.DAS2R;
phib = ( 84381.412819     + ...
       (   -46.811016     + ...
       (     0.0511268    + ...
       (     0.00053289   + ...
       (    -0.000000440  + ...
       (    -0.0000000176 ) ...
       * t) * t) * t) * t) * t) * Const.DAS2R;
psib = (    -0.041775     + ...
       (  5038.481484     + ...
       (     1.5584175    + ...
       (    -0.00018522   + ...
       (    -0.000026452  + ...
       (    -0.0000000148 ) ...
       * t) * t) * t) * t) * t) * Const.DAS2R;
epsa =  obl06(date1, date2);