function [FE_eq,FE_int,FE_fin,FE_prop] = FreeEnergyCmp(feql,fint,f,fprop,Nc,bc)

% Compare Free energies

% Equilbrium
% keyboard
[sigma,rho] = FEsigrhoCalcHR2D(feql,Nc);
FE_eq = FreeEnCalc(rho,sigma,bc);
% keyboard

% Initial
[sigma,rho] = FEsigrhoCalcHR2D(fint,Nc);
FE_int = FreeEnCalc(rho,sigma,bc);
% keyboard

% Final FE
[sigma,rho] = FEsigrhoCalcHR2D(f,Nc);
FE_fin = FreeEnCalc(rho,sigma,bc);
% keyboard

% FE from prop
[sigma,rho] = FEsigrhoCalcHR2D(fprop,Nc);
FE_prop = FreeEnCalc(rho,sigma,bc);
% keyboard

fprintf('FE_eq = %f\nFE_int = %f\nFE_fin = %f\nFE_prop = %f\n\n', FE_eq, FE_int, FE_fin, FE_prop);

end
