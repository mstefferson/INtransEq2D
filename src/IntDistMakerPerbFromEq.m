% IntDistMakerPerbFromEq.m
% Grid stuff
Nx    = 64;
Ny    = 64;
Nm    = 64;
L_box = 2*pi;
dx    = L_box / Nx;
dy    = L_box / Ny;
d_phi = 2*pi / Nm;

x     = 0  : dx : L_box-dx;
y     = x;
phi   = -pi:d_phi:pi-d_phi;
% Distribution stuff
Nc    = 10;            % Number of Coefficients
Modes = 3;
%Particle stuff
ParticleNum     = 100;
L_rod = 1;

b     = L_rod^2/pi;                       % Average excluded volume
c     = ParticleNum / (L_box);        % Concentration
bc    = b*c;

[Coeff_best, CoeffMat] = CoeffCalcExpCos2D(Nc,phi,bc);

f = DistBuilderExpCos2Dsing(Nc,phi,Coeff_best);

for ii = 1:Modes
    
    CoeffPert = 1 / (2*pi) ;
    f_perturbed = f + (CoeffPert + CoeffPert/10) + CoeffPert * cos(2*ii*phi); %Shift it so it's always positive
end

figure
subplot(2,1,1)
plot(phi,f)

subplot(2,1,2)
plot(phi,f_perturbed)

figure
subplot(2,1,1)
plot( real( fftshift( fft(f_perturbed)) ) )

subplot(2,1,2)
plot( imag( fftshift( fft(f_perturbed)) ) )



% Map this to a homogeneous system

rho = ones(Nx,Ny,Nm);
for i = 1:Nm
rho(:,:,i) = rho(:,:,i) .* f_perturbed(i);
end

CurrentNorm = trapz_periodic(y,trapz_periodic(x,trapz_periodic(phi,rho,3),2),1);
rho = rho .* ParticleNum ./ CurrentNorm;