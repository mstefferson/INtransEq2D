% Calculates the Nematic order parameter 
function [N, NOPx,NOPy] = NemOrderFnc( f, phi )
  
  cosPhi = cos(phi);
  sinPhi = sin(phi);

  Q_NOPxx_temp = trapz_periodic(phi, f .* (cosPhi .* cosPhi - 1/2) );
  Q_NOPxy_temp = trapz_periodic(phi, f .* (cosPhi .* sinPhi) );
  Q_NOPyy_temp = trapz_periodic(phi, f .* (sinPhi .* sinPhi - 1/2) );

  Q_temp = [Q_NOPxx_temp Q_NOPxy_temp; Q_NOPxy_temp Q_NOPyy_temp];
  [EigVec,EigS] = eigs(Q_temp);
  eigMaxQ_NOP = max(max(EigS));
  %Find the eigenvector corresponding to this eigenvalue
  NOPtemp = EigVec( EigS == repmat(max(EigS,[],2),[1,2]) );
  NOPx = NOPtemp(1);
  NOPy = NOPtemp(2);
  % We need to build 2x2 matrices and diagonalize them. But, we know how to
  % easily diagonalize a 2x2 matrix so can we find all the eigenvalues in 1
  % step.
  %  eigMaxQ_NOP = sqrt( ((Q_NOPxx_temp - Q_NOPyy_temp )./2).^2 + Q_NOPxy_temp.^2);
  % keyboard
  %Calculate Nematic Scalar parameter
  N = 2*eigMaxQ_NOP;

end

