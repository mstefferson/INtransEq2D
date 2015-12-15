% EqDistMakerExpLegMain.m
% Michael Stefferson
%
% Function returns the normalizated equilibrium distribution
% f(\theta) = exp( \sum a_2n P_2n (cos(\theta) ) / Z
% of hard rods. This is the distribution that minimizes the free energy.
%
% Inputs: 
% bc: scaled concentation (concent. * average excluded volume/particle)
% Nc: # of coefficients in the expansion
% Nx: Spatial resolution for integration
% plotme: turn on/off plotting of the distribution. Turn off if you're
% loopin!
%
% Output: 
%                   | a_2 0   0 ...  0   |
% CoeffTempMat =    | a_2 a_4 0 ...  0   |
%                   | ...           
%                   | a_2 a_4 ... a_{2Nc}| 
%
% f_mat =  [ f(a_2) ; f(a_2, a_4); ...; f(a_2, a_4,...a_{2Nc})

function [CoeffMat, f_best] = EqDistMakerMain2D(bc, Nc,Nx,plotme)
% bc  = 5;      % Scaled concentration. 
% Nc = 10;       % Number of coefficients.


% Integration stuff
% phi = linspace(-pi,pi,Nx);
phi = linspace(0,2*pi-2*pi/Nx,Nx);

% Calculate the coefficients of the expansion
[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,phi,bc);
% keyboard
% Build the distributions for various number of coefficients.
% f_mat    = DistBuilderExpCos2Dmat(Nc,phi,CoeffMat);
% Return the distribution made w/ the most coeff.
% f_best   = f_mat(Nc,:);  
f_best  = DistBuilderExpCos2Dsing(Nc,phi,Coeff_best);
% keyboard
% figure
% plot( theta, f_mat )
% title('Angular distribution')

% trapz(phi,f_best);
if plotme
figure
plot( phi, f_best )
titstr = sprintf('Angular distribution bc = %.2f', bc);
title(titstr)
ylabel('f(\phi)')
xlabel('\phi')
set(gca,'XLim',[min(phi) max(phi)] )
set(gca,'YLim',[0 max(f_best) + max(f_best)/ 10] )
end