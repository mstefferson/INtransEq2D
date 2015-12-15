% FEsigrhoCalcHR3D.m
%
% Calculates the rho and sigma in the free energy for 3D hard rods




function [sigma,rho] = FEsigrhoCalcHR3D(f,Nc)

Nt     = length(f);
x      = linspace(-1,1,Nt);
sigma  = 2*pi * trapz( x, f.*log(4.*pi.*f) );

d2nVec = KernCoeffCalcHardRod3D(Nc);

rho    = 1;
for n = 1:Nc
    Leg_2n_allm = legendre(2*n,x);
    Leg_2n      = Leg_2n_allm(1,:);
    rho        = rho - 16 * pi * d2nVec(n) * ( trapz(x,Leg_2n .* f) ) .^ 2; 
%     keyboard
end


end