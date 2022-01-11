function enu = EC2ENU(ece,orgece,orgllh)
%
%   EC2ENU: Convert ECEF coordinates to East-North-Up with
%   respect to orgece and orgllh (orgece is the same
%   location as orgllh).
%   Format of ece, orgece, and orgllh are: 1x3 matrix each

difece = ece - orgece;   % difference between coordinate
                         % origins

%   Rotate the difference vector into ENU coordinates

sla = sin(orgllh(1)); cla = cos(orgllh(1));
slo = sin(orgllh(2)); clo = cos(orgllh(2));

enu = [  -slo      clo      0 ; ...
       -sla*clo  -sla*slo  cla; ...
        cla*clo   cla*slo  sla] * difece';
enu = enu';