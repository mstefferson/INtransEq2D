% CoexSolveTest.m

function [Cn, Ci] = CoexSolveTest


alpha = 33.4;
rho   =  4 / ( sqrt( pi * alpha ) ) ;
sigma = log(alpha) - 1;
% If we want to pass a function in fsolve variables make a handle for it
% Then pass fsolve the handle

f_cn_handle = @(cn) NonLinFuncCn1(cn,rho,sigma);
Cn   =  fsolve(f_cn_handle,4);

f_ci_handle = @(ci) NonLinFuncCi1(ci,rho,sigma);
Ci    =  fsolve(f_ci_handle,4);

end

function F_ci = NonLinFuncCi1(ci,rho,sigma)
% Function we want to find the zero of in terms of ci
    F_ci = log( ( -1 + sqrt( 1 + 4 * rho * ci * ( 1 +  ci ) ) ) ./ ( 2*rho*ci ) ) ...
         - 1 + sqrt( 1 + 4 * rho * ci * ( 1 + ci ) ) + sigma - 2 * ci;
end

function F_cn = NonLinFuncCn1(cn,rho,sigma)
% Function we want to find the zero of in terms of cn
 F_cn = log( ( -1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) ) ./ ( 2*cn ) ) ...
         - 1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) - sigma - 2 * cn * rho;
end
