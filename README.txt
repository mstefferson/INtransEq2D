

% Functions 

[CoeffMat, f_best] = EqDistMakerMain3D(bc, Nc,Nx,plotme)
% Function returns the normalizated equilibrium distribution

[CoeffMat, f_best] = EqDistMakerMain2D(bc, Nc,Nx,plotme)
% Function returns the normalizated equilibrium distribution

PhaseTransFinder2D()
%Finds IN transition in 2D

PhaseTransFinder3D()
%Finds IN transition in 3D

In src:

[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,phi,bc)
% Calculates the coefficients for a 2D distribution of the form
% f(\theta) = exp( \sum a_2n cos(2*n*\theta) ) / Z

f_mat = DistBuilderExpCos2Dmat(Nc,phi,CoeffMat)
% Builds the normalized 3D distributions of the form
% f(\theta) = exp( \sum a_2n P_2n (cos(\theta) ) / Z
% f_mat =  [ f(a_2) ; f(a_2, a_4); ...; f(a_2, a_4,...a_{2Nc})

