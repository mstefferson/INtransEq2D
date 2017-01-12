% FreeEnCalc2D
% Michael Stefferson

function FE = FreeEnCalc(rho,sigma,bc)
FE = log(bc) + sigma + bc * rho;
end