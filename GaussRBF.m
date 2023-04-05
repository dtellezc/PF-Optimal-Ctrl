function Psi = GaussRBF( X,C, sig)
Psi = zeros(size(C, 2), size(X, 2));
Ctemp = C;
for i = 1:size(C, 2)
    C = repmat( Ctemp(:,i), 1, size(X,2) );
    psi = exp(-sum( (X - C).^2 ,1)/sig^2);
    Psi(i,:) = psi;
end
