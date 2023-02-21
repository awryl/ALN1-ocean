%Calcul de la matrice X'
[n,p1,p2] = size(F);
X = reshape(F,n,p1*p2);
Xbarre=mean(X);
Xprime = X - ones(n,1) * Xbarre;
Sigma = transpose(Xprime) * Xprime * (1/(nx-1));
fprintf('Computing EOFs...\n');
 %Les eofs sont dans la matrice EOF, et leur valeurs propres associés dans
 %la matrice diagonale D
 s=size(Sigma);
[ EOFs, lambda, droite ] = svd(transpose(Sigma));
fprintf('Done!\n');
fprintf('PCs...\n');
Pc = droite * lambda; % Principal components
fprintf('Done!\n');
fprintf('EigenValues!\n');
L = diag( lambda ) .^ 2 / (s(1)-1); % Valeurs propres
fprintf('Done!\n');