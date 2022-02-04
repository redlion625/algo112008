function gast = gst06a(uta,utb,tta,ttb)
%  Given:
%     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
%     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
%
%  Returned (function value):
%                double    Greenwich apparent sidereal time (radians)
%
%  Notes:
%
%  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
%     Julian Dates, apportioned in any convenient way between the
%     argument pairs.  For example, JD(UT1)=2450123.7 could be
%     expressed in any of these ways, among others:
%
%             uta            utb
%
%         2450123.7           0.0       (JD method)
%         2451545.0       -1421.3       (J2000 method)
%         2400000.5       50123.2       (MJD method)
%         2450123.5           0.2       (date & time method)
%
%     The JD method is the most natural and convenient to use in
%     cases where the loss of several decimal digits of resolution
%     is acceptable (in the case of UT;  the TT is not at all critical
%     in this respect).  The J2000 and MJD methods are good compromises
%     between resolution and convenience.  For UT, the date & time
%     method is best matched to the algorithm that is used by the Earth
%     rotation angle function, called internally:  maximum precision is
%     delivered when the uta argument is for 0hrs UT1 on the day in
%     question and the utb argument lies in the range 0 to 1, or vice
%     versa.

rnpb = pnm06a(tta,ttb);
gast = gst06(uta,utb,tta,ttb,rnpb);
