% CoexFinderINtrans2D.m


Nt = 5000;
phi  = linspace(-pi,pi,Nt);
Nc = 6;
epsilon = 10^(-5);

% Interation variables
CoeffVec      = rand(1,Nc);
CoeffVec_next = rand(1,Nc);
ci      = 4;
ci_next = 5;
cn      = 4;
cn_next = 5;

% Calculate the Legendre coefficients of the kernal
d2nVec          = KernCoeffCalcHardRod2D(Nc);
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
    f = DistBuilderExpCos2Dsing(Nc,phi,CoeffVec);
    [sigma,rho] = FEsigrhoCalcHR2D(f,Nc);
   
    %Calculate the next coefficients
    [cn_next, ci_next] = CoexConcCalcExpnsn(rho,sigma);
    
    CosSum = zeros(1, length(phi) );
    % Sum over all the Legendre polynomials
    for n = 1 : Nc
        CosSum   =  CosSum + CoeffVec(n) *  cos(2*n*phi);
    end % End of Legendre poly sums
    
    %     keyboard
    for j = 1:Nc
%         keyboard
        CoeffVec_next(j)  = pi * cn * d2nVec(j) ...
                * ( trapz( phi, cos(2*j*phi) .* exp(CosSum) ) ./ ...
                (trapz( phi, exp(CosSum) ) ) );
    end
    
    InterVarVec_next(1:Nc) = CoeffVec_next;
    InterVarVec_next(Nc+1) = ci_next;
    InterVarVec_next(Nc+2) = cn_next;
    if mod(counter, 10000) == 0
       keyboard 
    end
    counter = counter + 1;
%     keyboard
end % End coefficient interation





