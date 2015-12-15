% GaussAnsatzTester2D.m

% Spatial stuff
Nm    = 1000;
phi   = linspace(-pi,pi,Nm);

% Parameters
alpha = 17.93;

% Gaussians
G_old = sqrt(alpha/(2*pi)) * exp( - alpha * phi.^2 / 2);

figure
plot(phi,G_old)

trapz(phi,G_old)

% Split up phi

phi1 = phi(1:Nm/4);
phi2 = phi(Nm/4+1:3*Nm/4);
phi3 = phi(3*Nm/4+1:end);

G1   = 1/2*sqrt(alpha/(2*pi)) * exp( - alpha * (phi1+pi).^2 / 2);
G2   = 1/2*sqrt(alpha/(2*pi)) * exp( - alpha * (phi2).^2 / 2);
G3   = 1/2*sqrt(alpha/(2*pi)) * exp( - alpha * (phi3-pi).^2 / 2);

G_new = [G1 G2 G3];


figure
plot(phi,G_new)
trapz(phi,G_new)