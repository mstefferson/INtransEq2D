function NemPlotter( M, nxVec, bcVec, depParam )

numBc = length( bcVec );
numNx = length( nxVec );

if numNx > 1
  legCell = cell( 1, numBc );
  figure()
  hold
  for ii = 1:numBc
    plot( nxVec, M(:,ii) )
   legCell{ii} =  ['bc = '  num2str( bcVec(ii) ) ];
  end
  titStr = [ depParam ' vs $$N_{gr}$$'];
  title(titStr)
  xlabel('$$N_{gr}$$');
  ylabel(depParam);
  legH = legend( legCell );
  legH.Interpreter = 'latex';
  legH.Location = 'best';
end

if numBc > 1
  legCell = cell( 1, numNx );
  figure()
  hold
  for ii = 1:numNx
    plot( bcVec, M(ii,:) )
    legCell{ii} = ['$$N_{gr}$$ = '  num2str( nxVec(ii) ) ];
  end
  titStr = [depParam ' vs bc'];
  title(titStr)
  xlabel('bc');
  ylabel(depParam);
  legH = legend( legCell );
  legH.Interpreter = 'latex';
  legH.Location = 'best';
end

