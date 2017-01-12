% KernCoeffCalcHardRod3D.m
% Michael Stefferson
%
% Calculate the coefficients of the Legendre expansion of the kernal
% |sin \gamma| for a 3d hard rod interactions.

function d2nVec = KernCoeffCalcHardRod3D(Nc)
% Build a vector of the kernal's coefficients.
d2nVec        = zeros(1,Nc);

for n = 1:Nc;
    d2nVec(n) = ( pi * ( 4*n + 1 ) * doubfact( 2*n - 3 ) * doubfact( 2*n - 1) ) ...
           ./ ( 2 ^( 2*n + 2 ) * factorial(n) * factorial(n+1) );
end

% keyboard
end