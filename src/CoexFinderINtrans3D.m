% CoexFinderINtrans3D.m
addpath('C:\Users\MWS\Documents\MATLAB\Programs\MathRoutines\')

Nt = 5000;
x  = linspace(-1,1,Nt);
Nc = 9;
epsilon = 10^(-8);

% Interation variables
CoeffVec      = zeros(1,Nc);
CoeffVec_next = ones(1,Nc);
ci      = 4;
ci_next = 5;
cn      = 4;
cn_next = 5;


% Matrices
Leg_2n_Mat = zeros(Nc,Nt);
for i = 1:Nc
    Leg_2n_allm      =  legendre(2*i,x);
    Leg_2n           =  Leg_2n_allm(1,:);
    Leg_2n_Mat(i,:)  =  Leg_2n;
end
% Calculate the Legendre coefficients of the kernal
d2nVec   = KernCoeffCalcHardRod3D(Nc);
InterVarVec      = zeros(1,Nc+2);       %[a_2, a_4,... a_2Nc ci cn]
InterVarVec_next = ones(1,Nc+2);

% Iterate the coefficients
counter = 1;
while max(abs(InterVarVec - InterVarVec_next)) > epsilon
    
    % Update everything
    CoeffVec     = CoeffVec_next;
    ci           = ci_next;
    cn           = cn_next;
    InterVarVec  = InterVarVec_next;
    
    % Calculate the distribution, rho, and sigma at this step
    f = DistBuilderExpLeg3Dsing(Nc,x,CoeffVec);
    [sigma,rho] = FEsigrhoCalcHR3D(f,Nc);
  
%     keyboard
  
    %Calculate the next coefficients
    [cn_next, ci_next] = CoexConcCalcExpnsn(rho,sigma);
    
    LegPolSum = zeros(1, length(x) );
    % Sum over all the Legendre polynomials
    for n = 1 : Nc
        LegPolSum   =  LegPolSum + CoeffVec(n) *  Leg_2n_Mat(n,:);
    end % End of Legendre poly sums
    
    %     keyboard
    for j = 1:Nc
        %             keyboard
        CoeffVec_next(j) = 8 * cn * d2nVec(j) / pi  ...
            * ( trapz(x, Leg_2n_Mat(j,:) .* exp(LegPolSum)  ) ./ ...
            (trapz( x, exp(LegPolSum) ) ) );
    end
    
    InterVarVec_next(1:Nc) = CoeffVec_next;
    InterVarVec_next(Nc+1) = ci_next;
    InterVarVec_next(Nc+2) = cn_next;
    counter = counter + 1;
%     keyboard
end % End coefficient interation





