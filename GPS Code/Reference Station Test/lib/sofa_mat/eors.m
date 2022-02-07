function eo = eors(rnpb,s)
% 
% x = rnpb(2,0);
% ax = x / (1.0 + rnpb(2,2));
% xs = 1.0 - ax * x;
% ys = -ax * rnpb(2,1);
% zs = -x;
% p = rnpb(0,0) * xs + rnpb(0,1) * ys + rnpb(0,2) * zs;
% q = rnpb(1,0) * xs + rnpb(1,1) * ys + rnpb(1,2) * zs;
% eo = ((p != 0) || (q != 0)) ? s - atan2(q, p) : s;

x = rnpb(3,1);
ax = x / (1.0 + rnpb(3,3));
xs = 1.0 - ax * x;
ys = -ax * rnpb(3,2);
zs = -x;
p = rnpb(1,1) * xs + rnpb(1,2) * ys + rnpb(1,3) * zs;
q = rnpb(2,1) * xs + rnpb(2,2) * ys + rnpb(2,3) * zs;
if ((p ~= 0) || (q ~= 0))
    eo = s - atan2(q, p);
else
    eo = s;
end