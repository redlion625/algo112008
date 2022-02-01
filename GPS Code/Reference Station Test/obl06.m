function epsa = obl06(date1,date2)
t = ((date1 - Const.DJ00) + date2) / Const.DJC;

% Mean obliquity. 
epsa = (84381.406     + ...
       (-46.836769    + ...
       ( -0.0001831   + ...
       (  0.00200340  + ...
       ( -0.000000576 + ...
       ( -0.0000000434) * t) * t) * t) * t) * t) * Const.DAS2R;
