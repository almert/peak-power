function [out] = lentzs(nu,z)
% Y = steeds(nu,z) implements Lentz's method for computing the continued
% fraction, which gives the ratio of two modified Bessel functions
% evaluated at z. The output is I_{nu}(z)/I_{nu-1}(z)
% Here, nu is the dimension, z is the input of the Bessel functions.

if abs(z) <= 1e-30 % This is to avoid z^{-1} causing overflows
    out = z;
    return;
end

eps = 1e-10;
tiny = 1e-30;
fn = tiny;
Cn = fn;
Dn = 0;

inc = 2/z;
bn = (2*nu-2)/z;
while true
    bn = bn+inc;
    Dn = bn+Dn;
    if abs(Dn)<tiny
        Dn = tiny;
    end
    Cn = bn+1/Cn;
    if abs(Cn)<tiny
        Cn = tiny;
    end
    Dn = 1/Dn;
    df = Cn*Dn;
    fn = fn*df;
    if abs(df-1)<eps
        break;
    end
end
out = fn;
end

