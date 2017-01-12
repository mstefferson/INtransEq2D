% Wrapper for NemOrderFnc
function [nMat,fMax,fMin] = NemOrderDistMaxWrap( nxVec, bcVec, nC )
addpath('./src')
numNx = length( nxVec );
numBc = length( bcVec );
nMat = zeros( numNx, numBc ); 
fMax = zeros( numNx, numBc ); 
fMin = zeros( numNx, numBc ); 

for ii = 1:numNx
  for jj = 1:numBc
    [~, f_best, phi] = EqDistMakerMain2D( bcVec(jj), nC, nxVec(ii), 0 );
    fMax(ii,jj) = max( f_best );
    fMin(ii,jj) = min( f_best );
    nMat(ii,jj) = NemOrderFnc( f_best, phi );
  end
end
  
end

