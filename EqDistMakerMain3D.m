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

function [CoeffMat, f_best] = EqDistMakerMain3D(bc, Nc,Nx,plotme)
addpath('./src')
% bc  = 5;      % Scaled concentration. 
% Nc = 10;       % Number of coefficients.


% Integration stuff
dx = 2 / Nx;
x = -1 : dx : 1 - dx;

% Calculate the coefficients of the expansion
% CoeffMat = CoeffCalcExpLeg3D_test(Nc,x,bc,d2nVec);
[Coeff_best, CoeffMat] = CoeffCalcExpLeg3D(Nc,x,bc);
% keyboard
% Build the distributions for various number of coefficients.
f_mat    = DistBuilderExpLeg3Dmat(Nc,x,CoeffMat);
% Return the distribution made w/ the most coeff.
f_best   = f_mat(Nc,:);  
% figure
% plot( theta, f_mat )
% title('Angular distribution')

if plotme
figure
theta = linspace(0,pi,Nx);
plot( theta, f_best )
titstr = sprintf('Angular distribution bc = %.2f', bc);
title(titstr)
ylabel('f(\theta)')
xlabel('\theta')
set(gca,'XLim',[0 pi] )
set(gca,'YLim',[0 max(f_best) + max(f_best)/ 10] )
end
