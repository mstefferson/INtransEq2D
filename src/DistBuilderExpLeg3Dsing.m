% DistBuilderExpLeg3D.m
% Michael Stefferson
%
% Builds the normalized 3D distributions of the form
% f(\theta) = exp( \sum a_2n P_2n (cos(\theta) ) / Z
%
% Inputs:
% Nc: Max number of coeff
% x : spatial coordin. vector. x = cos(\theta)
% CoeffMat: Matrix of coefficients
%
% Outputs:
% f_mat =  [ f(a_2) ; f(a_2, a_4); ...; f(a_2, a_4,...a_{2Nc})

function f = DistBuilderExpLeg3Dsing(Nc,x,Coeff)
% Distribution
f = zeros(1,length(x));
LegPolSum     = zeros(1, length(x) );
for i = 1:Nc  
          Leg_k_allm    =  legendre(2*i,x);
          Leg_k         =  Leg_k_allm(1,:);
          LegPolSum     = LegPolSum + Coeff(i) * Leg_k;    
%           keyboard
end
%        Normalize it again to be safe.
       f = exp( LegPolSum ) ./ (2*pi * trapz( x, exp( LegPolSum ) ) );
% keyboard