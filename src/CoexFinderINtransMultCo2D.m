% CoexFinderINtransMultCo2D.m
addpath('C:\Users\MWS\Documents\MATLAB\Programs\MathRoutines\')

Nx       = 5000;
x        = linspace(-1,1,Nx);
epsilon  = 10^(-5);
Nc_End   = 8;
Nc_Start = 3;
Nc_Diff  = Nc_End-Nc_Start + 1;


% Calculate the Legendre coefficients of the kernal

AllCoeffMat = zeros(Nc_Diff, Nc_End+3);
AllCoeffMat(:,1) = Nc_Start:Nc_End;
d2nVec   = KernCoeffCalcHardRod2D(Nc_End);

MatPlaceHolder = 1;
for Nc = Nc_Start:Nc_End
    % Interation variables
    %     keyboard
    CoeffVec      = zeros(1,Nc);
    CoeffVec_next = ones(1,Nc);
    ci      = 4;
    ci_next = 5;
    cn      = 4;
    cn_next = 5;
    
    
    InterVarVec      = zeros(1,Nc+2);       %[a_2, a_4,... a_2Nc ci cn]
    InterVarVec_next = ones(1,Nc+2);
    
    % Iterate the coefficients
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
%             keyboard
        end
        counter = counter + 1;
        %     keyboard
    end % End coefficient interation
    
    %     keyboard
    % Update Matrix
    
    AllCoeffMat(MatPlaceHolder,2:Nc+1)   = CoeffVec_next;
    AllCoeffMat(MatPlaceHolder,Nc_End+2) = ci_next;
    AllCoeffMat(MatPlaceHolder,Nc_End+3)   = cn_next;
    MatPlaceHolder = MatPlaceHolder + 1;
end



