% CoexSandBox.m

function [] = CoexSandBox
Ci_co1  = 3.45;
Cn_co1  = 5.12;
alpha_1 = 33.4;

Ci_co2  = 0.98;
Cn_co2  = 0.65;
alpha_2 = 0.538;


rho1   =  4 / ( sqrt( pi * alpha_1 ) ) ;
sigma1 = log(alpha_1) - 1;

Cn1_calc = ( -1 + sqrt( 1 + 4 * rho1 * Ci_co1 * ( 1 + Ci_co1 ) ) ) ./ (2*rho1);

Ci1_calc = ( -1 + sqrt( 1 + 4 * Cn_co1 * ( 1 + rho1 * Cn_co1) ) ) ./2;

rho2 =  4 / ( sqrt( pi * alpha_2 ) ) ;
sigma2 = log(alpha_2) - 1;

Cn2_calc = ( -1 + sqrt( 1 + 4 * rho2 * Ci_co2 * ( 1 + Ci_co2 ) ) ) ./ (2*rho2);
Ci2_calc = ( -1 + sqrt( 1 + 4 * Cn_co2 * ( 1 + rho2 * Cn_co2) ) )  ./2;

F_cn1  = NonLinFuncCn1(Cn_co1,rho1,sigma1);
disp(F_cn1)

%Use solve to find zeros of a function
cn_zero = fsolve(@NonLinFuncCn2,4);
f_cn_handle = @(cn) NonLinFuncCn1(cn,rho1,sigma1);
% Pass the handle to fsolve.
cn_zero2    =  fsolve(f_cn_handle,4);

% See how close to zero solution is
F_ci1       = NonLinFuncCi1(Ci_co1,rho1,sigma1);

% Find zero without passing it variables
ci_zero = fsolve(@NonLinFuncCi2,4);

% If we want to pass the function variables make a handle for it
f_ci_handle = @(ci) NonLinFuncCi1(ci,rho1,sigma1);
% Pass the handle to fsolve.
ci_zero2    =  fsolve(f_ci_handle,4);

keyboard
F_cn2 = NonLinFuncCn1(Cn_co2,rho2,sigma2);
disp(F_cn2)

F_ci2 = NonLinFuncCi1(Ci_co2,rho2,sigma2);
disp(F_ci2)

%
Ci = ci_zero;
Cn = cn_zero;

% Coexistence equations
CoEx1LHS = Ci*(1+Ci);
CoEx1RHS = Cn*(1+rho*Cn);

CoEx2LHS = log(Ci) + 2*Ci;
CoEx2RHS = log(Cn) + sigma + 2*Cn*rho;

keyboard
end

function F_cn = NonLinFuncCn1(cn,rho,sigma)
% Function we want to find the zero of in terms of cn
 F_cn = log( ( -1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) ) ./ ( 2*cn ) ) ...
         - 1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) - sigma - 2 * cn * rho;
end

function F_cn = NonLinFuncCn2(cn)

alpha = 33.4;
% alpha = 0.538;
rho   =  4 / ( sqrt( pi * alpha ) ) ;
sigma = log(alpha) - 1;

% Function we want to find the zero of in terms of cn
 F_cn = log( ( -1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) ) ./ ( 2*cn ) ) ...
         - 1 + sqrt( 1 + 4 * cn * ( 1 + rho * cn ) ) - sigma - 2 * cn * rho;
end


function F_ci = NonLinFuncCi1(ci,rho,sigma)
% Function we want to find the zero of in terms of ci
    F_ci = log( ( -1 + sqrt( 1 + 4 * rho * ci * ( 1 +  ci ) ) ) ./ ( 2*rho*ci ) ) ...
         - 1 + sqrt( 1 + 4 * rho * ci * ( 1 + ci ) ) + sigma - 2 * ci;
end


function F_ci = NonLinFuncCi2(ci)

% alpha = 33.4;
% alpha = 0.538;
rho   =  4 / ( sqrt( pi * alpha ) ) ;
sigma = log(alpha) - 1;

% Function we want to find the zero of in terms of ci
    F_ci = log( ( -1 + sqrt( 1 + 4 * rho * ci * ( 1 +  ci ) ) ) ./ ( 2*rho*ci ) ) ...
         - 1 + sqrt( 1 + 4 * rho * ci * ( 1 + ci ) ) + sigma - 2 * ci;
end

