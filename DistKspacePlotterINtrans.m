% DistKspacePlotterINtrans.m

bc = 4.5;
Dim = 2;
Nx = 64;
Nc = 10;

if Dim == 2
  %Length of interval  
   L = 2*pi;
 %Build the distribution
   [CoeffMat, f] = EqDistMakerMain2D(bc, Nc,Nx,0); 
   x = linspace(-pi,pi,Nx);
else 
   L = pi;  
  [CoeffMat, f] = EqDistMakerMain3D(bc, Nc, Nx, 0);
  x = linspace(0,pi,Nx);
end

 % k-space stuff  
   dx = L/Nx;
   kx_max   = pi / dx;
   delta_kx = 2*pi / L;
   kx = [-kx_max: delta_kx: kx_max-delta_kx];
%    keyboard
f_FT = fftshift(fft(f));

figure
subplot(3,1,1)
plot( x,f )
subplot(3,1,2)
plot( kx,real( f_FT ))
subplot(3,1,3)
plot( kx,imag( f_FT ))

