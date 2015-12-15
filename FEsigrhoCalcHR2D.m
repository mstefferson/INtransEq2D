

function [sigma,rho] = FEsigrhoCalcHR2D(f,Nc)
% Turn f into a rho vector
[n, m]  = size(f);
if n > m
   f = f'; 
end
Nt     = length(f);
phi    = linspace(-pi,pi,Nt);
sigma  = trapz( phi, f.*log(2.*pi.*f) );

d2nVec = KernCoeffCalcHardRod2D(Nc);

% Calculate rho using the expansion of the kernal
rho    = 1;
% keyboard
for n = 1:Nc
%     keyboard
    rho = rho -  pi * d2nVec(n) / 2 * ( trapz( phi, cos(2*n*phi) .* f ) ) .^ 2; 
end

if 0
% Calculate rho doing the double integral
DistProd =  f' * f;
[phi2, phiprime2] = meshgrid(phi,phi);

PhiDiff = phi2 - phiprime2; 
Kernal    =  abs( sin(PhiDiff) );
Integrand =  DistProd .* Kernal;

rhoDoubInt = pi / 2 * trapz( phi, trapz(phi,Integrand,1), 2);
end
% keyboard
end