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

function CoeffMat = CoeffCalcExpLeg3D_test(Nc,x,bc,d2nVec)
 
% Calculate the coefficents of the distribution
CoeffMat = zeros(Nc);
% Condition for while loop. How much coefficient need to be changing by 
% before stopping
epsilon = 1e-12;        
for i = 1:Nc
   
    CoeffTempVec_next = CoeffMat(i,:);
    CoeffTempVec      = ones(1,length(CoeffMat(i,:)));
    
    % Iterate the coefficients
    while max(abs(CoeffTempVec - CoeffTempVec_next)) > epsilon
        CoeffTempVec = CoeffTempVec_next;
        LegPolSum = zeros(1, length(x) );
        % Sum over all the Legendre polynomials
        for n = 1 : i
            Leg_2n_allm =  legendre(2*n,x);
            Leg_2n      =  Leg_2n_allm(1,:);
            LegPolSum  =  LegPolSum + CoeffTempVec(n) * Leg_2n;
            if i == 2
                %            keyboard
            end
        end % End of Legendre poly sums

        CoeffTempVec_next(1:i) = 4 * bc / pi * d2nVec(1:i) ...
            * ( 3 * trapz(x, x.^2 .* exp(LegPolSum) ./ ...
                 trapz(x, exp(LegPolSum) ) ) - 1 );

    end % End coefficient interation
    
%     keyboard
    % Update the matrix
    CoeffMat(i,:)   = CoeffTempVec_next;
    
    % Make a better guess
    if i ~= Nc
        CoeffMat(i+1,:) = CoeffTempVec_next;
    end
    
end % End loop over all coefficients.

% keyboard
end