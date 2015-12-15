% CoeffCalcExpLeg3D.m
% Michael Stefferson
%
% Calculates the coefficients for a distribution of the form
% f(\theta) = exp( \sum a_2n P_2n (cos(\theta) ) / Z
% Inputs:
% Nc: Number of max coefficients
% x : spatial coordin. vector. x = cos(\theta)
%
% Outputs:
%                   | a_2 0   0 ...  0   |
% CoeffTempMat =    | a_2 a_4 0 ...  0   |
%                   | ...            0
%                   | a_2 a_4 ... a_{2Nc}|

function [Coeff_best, CoeffMat] = CoeffCalcExpLeg3D(Nc,x,bc)


% Calculate the Legendre coefficients of the kernal
d2nVec   = KernCoeffCalcHardRod3D(Nc);
% Calculate the coefficents of the distribution
CoeffMat   = zeros(Nc);
Leg_2n_Mat = zeros(Nc, length(x));
% Condition for while loop. How much coefficient need to be changing by
% before stopping
epsilon = 1e-12;

for i = 1:Nc
    % Make initial guess for the coeff the previous coeff.
    CoeffTempVec_next = CoeffMat(i,:);
    CoeffTempVec      = ones(1,length(CoeffMat(i,:)));
    
    Leg_2n_allm      =  legendre(2*i,x);
    Leg_2n           =  Leg_2n_allm(1,:);
    Leg_2n_Mat(i,:)  =  Leg_2n;
    
    % Iterate the coefficients
    while max(abs(CoeffTempVec - CoeffTempVec_next)) > epsilon
        CoeffTempVec = CoeffTempVec_next;
        LegPolSum = zeros(1, length(x) );
        % Sum over all the Legendre polynomials
        for n = 1 : i
            LegPolSum   =  LegPolSum + CoeffTempVec(n) *  Leg_2n_Mat(n,:);
        end % End of Legendre poly sums
        
        for j = 1:i
%             keyboard
            CoeffTempVec_next(j) = 8 / pi * bc * d2nVec(j) ...
                * ( trapz(x, Leg_2n_Mat(j,:) .* exp(LegPolSum)  ) ./ ...
                  (trapz( x, exp(LegPolSum) ) ) );
        end
%         keyboard
    end % End coefficient interation
    
%     keyboard
    % Update the matrix
    CoeffMat(i,:)   = CoeffTempVec_next;
    
    % Make a better guess
    if i ~= Nc
        CoeffMat(i+1,:) = CoeffTempVec_next;
    end
    
end % End loop over all coefficients.
Coeff_best = CoeffMat(Nc,:);
% keyboard
end