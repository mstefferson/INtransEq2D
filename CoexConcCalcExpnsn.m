function [Cn, Ci] = CoexConcCalcExpnsn(rho,sigma)

% If we want to pass a function in fsolve variables make a handle for it
% Then pass fsolve the handle
%Suppress display message
options = optimset('Display','off');
f_cn_handle = @(cn) NonLinFuncCnExp(cn,rho,sigma);
Cn1   =  fsolve(f_cn_handle,4,options);
Ci1 = ( -1 + sqrt( 1 + 4 * Cn1 * ( 1 + rho * Cn1) ) ) ./2;


% f_ci_handle = @(ci) NonLinFuncCiExp(ci,rho,sigma,options);
% Ci2    =  fsolve(f_ci_handle,4);
% Cn2 = ( -1 + sqrt( 1 + 4 * rho * Ci2 * ( 1 + Ci2 ) ) ) ./ (2*rho);
% 
% keyboard
Cn = Cn1;
Ci = Ci1;
% 
% 
% % Coexistence equations
% CoEx1LHS = Ci*(1+Ci);
% CoEx1RHS = Cn*(1+rho*Cn);
% 
% CoEx2LHS = log(Ci) + 2*Ci;
% CoEx2RHS = log(Cn) + sigma + 2*Cn*rho;

% keyboard
end

function F_ci = NonLinFuncCiExp(ci,rho,sigma)
% Function we want to find the zero of in terms of ci
    F_ci = log( ( -1 + sqrt( 1 + 4 * rho * ci * ( 1 +  ci ) ) ) ./ ( 2*rho*ci ) ) ...
         - 1 + sqrt( 1 + 4 * rho * ci * ( 1 + ci ) ) + sigma - 2 * ci;
end

function F_cn = NonLinFuncCnExp(cn,rho,sigma)
% Function we want to find the zero of in terms of cn
 F_cn = log( ( -1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) ) ./ ( 2*cn ) ) ...
         - 1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) - sigma - 2 * cn * rho;
end
