% PhaseTransFinder2D
clear
% Add paths just in case
addpath('C:/Users/MWS/Documents/MATLAB/Programs/MathRoutines/');

% Grid, vec, and coefficient things
trial = 2;
N  = 10;          % How many different concentrations to use
Nx = 1000;        % Spatial grid points
Nc = 15;          % # Coefficients
CoeffDiff = 0.01;  % The amount the coefficients need to change to consider it a phase change

% Movie things
nFrames = N;
FPS     = 2;
FileNameMov = sprintf('PhaseTrans2D_t%d.avi',trial);
FileNameRhoSig = sprintf('RhoSigVsCon2D_t%d',trial);
FileNameDist = sprintf('DistPhaseTrans2D_t%d',trial);
x = linspace(-pi,pi,Nx);

% Concentration vector
bc_start = 1;
bc_end   = 3;
bc = linspace(bc_start,bc_end,N);

% Initialize some things
Coeff      = zeros(1,Nc);
Coeff_new  = zeros(1,Nc);
f          = zeros(1,Nx);
f_new      = zeros(1,Nx);
f_rec      = zeros(Nc,Nx);
rhoVec = zeros(1,N);
sigVec = zeros(1,N);
GotPhaseTrans = 0;

for i = 1:length(bc)
    % Calculate the distribution, rho, and sigma at a specific concentration
    [Coeff_new, f_new]     = EqDistMakerMain2D(bc(i), Nc, Nx, 0);
    [sigVec(i), rhoVec(i)] = FEsigrhoCalcHR2D(f_new,Nc);
    
    % See when the coefficients change. Skip the first one
    if GotPhaseTrans ~= 1
        if i ~= 1
            if abs( Coeff(Nc,1)-Coeff_new(Nc,1) ) > CoeffDiff
                %             keyboard
                bc_iso = bc(i-1);
                bc_nem = bc(i);
                
                f_iso = f;
                f_nem = f_new;
                GotPhaseTrans = 1;
            end
        end
    end
    %Update
    f     = f_new;
    Coeff = Coeff_new;
    %Store
    f_rec(i,:) = f_new;
end

if GotPhaseTrans
    fprintf('Phase transition between c'' = %f and c'' = %f \n', bc_iso,bc_nem);
else
    fprintf('Sorry my lord. I did not find a phase transition. \n')
    fprintf('Possibly consider making the phase transition condition a little more lenient\n')
    bc_iso = 0;
    bc_nem = 0;
    f_iso = 0;
    f_nem = 0;
end
% Plot everything
PhaseTransPlotter(bc,sigVec,rhoVec,x, bc_iso, bc_nem,f_iso,f_nem,...
    f_rec,nFrames,FileNameMov,FileNameDist, FileNameRhoSig,FPS,GotPhaseTrans)