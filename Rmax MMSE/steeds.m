function [out] = steeds(nu,z)
% Y = steeds(nu,z) implements Steed's method for computing the continued
% fraction, which gives the ratio of two modified Bessel functions
% evaluated at z. The output is I_{nu}(z)/I_{nu-1}(z)
% Here, nu is the dimension, z is the input of the Bessel functions.

if abs(z) <= 1e-30 % This is to avoid z^{-1} causing overflows
    out = z;
    return;
end

eps = 1e-10;
Dn = z/(2*nu);
dCn = Dn;
Cn = Dn;

inc = 2/z;
bn = 2*nu/z;
while true
    bn = bn+inc;
    Dn = 1/(Dn+bn);
    dCn = (bn*Dn-1)*dCn;
    Cn = Cn+dCn;
    if abs(dCn) < eps
        break;
    end
end
out = Cn;
end

