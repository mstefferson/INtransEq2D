function [Ca, Ci,alpha] = CoexConcCalcGauss(dim)

% If we want to pass a function in fsolve variables make a handle for it
% Then pass fsolve the handle
%Suppress display message
options = optimset('Display','off');

if dim == 2;
f_ci_handle = @(ci) NonLinFuncCi1Gauss2d(ci);
Ci          =  fsolve(f_ci_handle,4,options);
Ca          = Ci * ( 1 + Ci ) / 2;
alpha       = pi * Ca^2;
end

if dim == 3;
    
f_ci_handle = @(ci) NonLinFuncCi1Gauss3d(ci);
Ci          =  fsolve(f_ci_handle,4,options);
Ca          =  Ci * ( 1 + Ci ) / 3;
alpha       =  4 * Ca^2 / pi;
    
end
% 
% % Coexistence equations
% CoEx1LHS = Ci*(1+Ci);
% CoEx1RHS = Cn*(1+rho*Cn);
% 
% CoEx2LHS = log(Ci) + 2*Ci;
% CoEx2RHS = log(Cn) + sigma + 2*Cn*rho;

% keyboard
end

function F_ci = NonLinFuncCi1Gauss2d(ci)
% Function we want to find the zero of in terms of ci
    F_ci = 2 * ci - log( ci * ( 1 + ci ) ^ 2) -3/2 - log ( pi / (4 * sqrt(2)) );
end

function F_ci = NonLinFuncCi1Gauss3d(ci)
% Function we want to find the zero of in terms of ci
    F_ci = 2 * ci - log( ci^2 * ( 1 + ci ) ^ 3) -3 - log ( 4 / (27*pi) );
end