% PhaseTransFinder3D
clear
% Add paths just in case
addpath('C:/Users/MWS/Documents/MATLAB/Programs/MathRoutines/');

% Grid, vec, and coefficient things
trial = 3;
N  = 20;          % How many different concentrations to use
Nx = 10000;        % Spatial grid points 
Nc = 15;          % # Coefficients
CoeffDiff = 0.1;  % The amount the coefficients need to change to consider it a phase change


% Movie things
nFrames = N;
FPS     = 2;
FileNameMov = sprintf('PhaseTrans3D_t%d.avi',trial);
FileNameRhoSig = sprintf('RhoSigVsCon3D_t%d',trial);
FileNameDist = sprintf('DistPhaseTrans3D_t%d',trial);
x = linspace(-1,1,Nx);

% Concentration vector
bc_start = 2;
bc_end   = 5;
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
    [Coeff_new, f_new]     = EqDistMakerMain3D(bc(i), Nc, Nx, 0);
    [sigVec(i), rhoVec(i)] = FEsigrhoCalcHR3D(f_new,Nc);
    
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
end % End loop over concentrations

if GotPhaseTrans
    fprintf('Phase transition between c'' = %f and c'' = %f \n', bc_iso,bc_nem);
else
    fprintf('Sorry my lord. I did not find a phase transition.')
    fprintf('Possibly consider making the phase transition condition a little more lenient')
end

% Plot everything
PhaseTransPlotter(bc,sigVec,rhoVec,x, bc_iso, bc_nem,f_iso,f_nem,...
    f_rec,nFrames,FileNameMov,FileNameDist, FileNameRhoSig,FPS,GotPhaseTrans)
