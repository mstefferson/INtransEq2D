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

function f_mat = DistBuilderExpCos2Dmat(Nc,phi,CoeffMat)
% Distribution
f_mat = zeros(Nc,length(phi));
for i = 1:Nc   
       CosSum = zeros(1, length(phi) );
       for j = 1 : i
           CosSum     = CosSum + CoeffMat(i,j) * cos(2*j*phi);  
       end
%        Normalize it again to be safe.
       f_mat(i,:) = exp( CosSum ) ./ ( trapz( phi, exp( CosSum ) ) );
end

% keyboard