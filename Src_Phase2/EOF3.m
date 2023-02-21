function EOF(datafile,doplot,threshold,algo)
% datafile = name of the file containing data or 'gendata' to use the synthetic data generator
% doplot = true if plots are requested
% threshold = threshold for filtering eigenvalues, in terms of percentage of variability to be recovered
% algo = 0,1,2 : versions of the fortran algorithms / 3 : matlab implementation

% in .cshrc and .login files :
% source /mnt/n7fs/ens/tp_guivarch/opt/intel/Compiler/11.0/083/bin/ifortvars.csh intel64
% setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH":/mnt/n7fs/ens/tp_guivarch/opt/libstd/usr/lib/"
% setenv LD_RUN_PATH $LD_LIBRARY_PATH

close all;

!ifort -c -O -fPIC -vec-report0 m_subspace_iter.f90
mex mex_subspace_iter.F90 subspace_iter.f90 m_subspace_iter.f90  FC="ifort" FFLAGS="-O -fPIC -vec-report0" CC="" LD="ifort" FLIBS="-L/applications/matlabr2011a/bin/glnxa64 -lmx -lmex -lmat -lm -L/mnt/n7fs/ens/tp_guivarch/opt/OpenBLAS/ -lopenblas";

if (algo<0 | algo>3)
  algo = 3
end

fprintf('......................................\n');
fprintf('  LOADING/GENERATING DATA \n');
fprintf('......................................\n\n');

if(~strcmp(datafile,'gendata'))
 % Load data file; the field is in F, coordinates are in lat and lon
 fprintf('Loading data...\n');
 load(datafile);
else
 % Use the synthetic data generator
 nx=24;
 ny=24;
 nt=200;
 neof=4;
 addnoise=false;
 mydata=gendata(nx,ny,nt,neof,addnoise);
 lon=1:nx;
 lat=1:ny;
 time=1:nt;
end

nx=length(lon);
ny=length(lat);
nt=length(time);

%% 1. EOF on the whole data minus last year.

% Apres avoir trace les PCs on constate une grosse anomalie sur le debut des donnees de Kaplan, on tronque donc les 80 premiers mois
if (strcmp(datafile,'Kaplan.mat'))
 F=mydata(80:nt-12,:);
else
 F=mydata(1:nt-12,:);
end

% Remove undefined values
F(F==-32768)=0;
% Animated plot (can be very long, use with care...)
startplot=input('Skip plotting data ? (ENTER : Yes , 0 : No)');
if(doplot & startplot==0)
  figure; 
  fprintf('ENTER to go on. \n');
  pause;
  fprintf('Displaying... \n');
  for i=1:nt-12
    mi=min(min(F(i,:)));
    ma=max(max(F(i,:)));
    if(mi~=ma)
      myfig = reshape(F(i,:),ny,nx);
      %contour(lon,lat,myfig,mi:(ma-mi)/30:ma);
      imagesc(myfig(end:-1:1,:),[mi ma]);
      drawnow;
    end
    if(mod(i,100)==0)
      stopplot=input('0 to stop, ENTER to go on.');
      if(stopplot==0)
        break;
      end
    end
  end
end


fprintf('......................................\n');
fprintf('          DATA ANALYSIS       \n');
fprintf('......................................\n');

%% 2. Calcul de Z
fprintf('\n Computing anomaly matrix and covariance matrix...\n');
Z = detrend(F);

%% 3. Calcul de S
S = Z'*Z;
fprintf('\n Done.\n');

%% 4. Recherche des paires propres de S



% Algorithme fortran v0
if algo == 0
  fprintf('\n Computing EOFs... (V0 Fortran)\n');
  v0=rand(nx*ny,12);
  maxit = 5000;
  epsi=1e-12;
  tic;
  [V,vp,res_ev,it_ev,it,nd] = mex_subspace_iter(0,S,1,v0,threshold/100,maxit,epsi);
  toc;
end
     
% Algorithme fortran v1
if algo == 1
  fprintf('\n Computing EOFs... (V1 Fortran)\n');
  v0=rand(nx*ny,12);
  maxit = 5000;
  epsi=1e-12;
  tic;
  [V,vp,res_ev,it_ev,it,nd] = mex_subspace_iter(1,S,1,v0,threshold/100,maxit,epsi);
  toc;
end
     
% Algorithme fortran v2
if algo == 2
  fprintf('\n Computing EOFs... (V2 Fortran)\n');
  v0=rand(nx*ny,12);
  maxit = 5000;
  epsi=1e-12;
  tic;
  [V,vp,res_ev,it_ev,it,nd] = mex_subspace_iter(2,S,3,v0,threshold/100,maxit,epsi);
  toc;
end

% Algorithme matlab

if algo == 3

  fprintf('\n Computing eigenpairs...\n');
  tic;
  [EOF,D]=eig(S);
  toc;
  fprintf('\n Done.\n');

  % Determination du nombre de valeurs propres à retenir nd
  % On retient var = PercentTrace% de la variance

  fprintf('\n Keeping important eigenpairs and displaying their dominance...\n');

  % Pour eviter des calculs trop longs voire infinis si on donne threshold > 100
  if threshold > 99.9999
    threshold = 99.9999
  end

  PercentTrace=threshold/100;
  nd=1;
  TraceDeD = trace(D);
  % La plus grande valeur propre est a la fin de D. On aurait pu la "renverser" mais apres test cela prend presque autant de temps que le calcul des eigenpairs..
  p=size(D,1);
  dominance = D(p+1-nd,p+1-nd)/TraceDeD
  TotalPercent = dominance;
  while TotalPercent < PercentTrace
    % Calcul de la dominance de la valeur propre 
    nd = nd+1;
    dominance = D(p+1-nd,p+1-nd)/TraceDeD
    TotalPercent = TotalPercent + dominance;
  end

  % On ne retient que les nd premiers EOF, qu'on met dans V (comme pour D, on commence par la fin)
  for i=1:nd
    V(:,i)=EOF(:,p+1-i);
  end

end

fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

%% 5. Tracer les EOF retenus et leur PC

fprintf('\n Plotting EOFs & PCs...\n');

% Plot des EOF
if(doplot)
  figure;
  fprintf('Displaying...\n');
  for i=1:nd
    % On veut plot les EOF à gauche, les PC au milieu et un zoom sur les PC a droite. Pour les EOF, il faut donc que leur position sur le subplot soit congru à 1[3]
    subplot(nd,3,i*3-2);
    mi=min(min(V(:,i)));
    ma=max(max(V(:,i)));
    if(mi~=ma)
      myfig = reshape(V(:,i),ny,nx);
      %contour(lon,lat,myfig,mi:(ma-mi)/30:ma);
      imagesc(myfig(end:-1:1,:),[mi ma]);
    end
    if(mod(i,100)==0)
      stopplot=input('0 to stop, ENTER to go on.\n');
      if(stopplot==0)
        break;
      end
    end
  end
end

% Plot des PC
PC = Z*V;
if(doplot)
  fprintf('Displaying...\n');
  for i=1:nd
    % On veut les PC au milieu, avec un zoom a droite (sinon la courbe est inexploitable)
    subplot(nd,3,i*3-1);
    plot(PC(:,i));
    subplot(nd,3,i*3);
    % Pour zoomer, on prend quelques donnees situees vers le milieu de la plage de donnees
    plot(PC(floor(nt/2):floor(nt/2 + nt/50),i)); 
  end
  drawnow;
end

fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

%% 6. Verifier la qualite de la base : ||Z.V.V' - Z||/||Z|| doit etre petit
fprintf('\n Computing the quality of the EOF basis (the error should be low) :\n\n');
errorofv = norm(Z*V*V' - Z)/norm(Z)
fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

%% 7. Verifier que la qualite est meilleure que celle d'une base aleatoire
fprintf('\n    Comparing the previously computed quality with the quality of a random basis...');
fprintf('\n Generating a random basis... \n');
ns = nx*ny;
W = rand(ns,nd);
% On continue a generer des vecteurs jusqu'a ce que ca soit effectivement une base
while rank(W)<nd
  W = rand(ns,nd);
end
fprintf('Done.\n');
fprintf('\n Computing the quality of the random basis (the error should be higher) :\n\n');
errorrandom = norm(Z*W*W' - Z)/norm(Z)
fprintf('\n Done.\n');

fprintf('\n......................................\n');
fprintf('       DATA PREDICTION      \n');
fprintf('........................................\n');

% On a deja calcule les EOF de F sans la derniere annee dans la partie analyse

fprintf('\n Shaping data for calculation... \n');
Fy = mydata(nt-11:nt,:);
Fy(Fy==-32768)=0;

% On choisit (arbitrairement) une partie connue
% exemple : knowndata = floor(ns/4):floor(3*ns/4);
knowndata = 1:(ns-1);

% On restreint Fy et V à cette partie connue (on appelle ces restrictions knownFy et Vi)
knownFy = mydata(nt-11:nt,knowndata);
knownFy(knownFy==-32768)=0;

Vi = V(knowndata,:);

% Si les deux sont differents, la decomposition QR ne pourra pas se faire convenablement
nd
rankofVi = rank(Vi)

fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

fprintf('\n QR decomposition of the restricted EOFs... \n');
[Q,R] = qr(Vi);
fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

fprintf('\n Solving the least-square problem... \n');
% On doit resoudre Vi*alpha = Fy' avec Vi = Q*R => Q*R*alpha = Fy' => alpha = inv(R) * Q' * Fy'
alpha = R\(Q'*knownFy');
fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

Fp=(V*alpha)';

fprintf('\n Computing the quality of the prediction (the error should be low)... \n');
% Calcul de la qualite de la prediction :
errorpredic = norm(Fp-Fy)/norm(Fy)
fprintf('\n Done.\n');

fprintf('\nENTER to go on. \n');
pause;

% On refait les calculs avec une base aleatoire (on utilise la meme que pour l'analyse)
fprintf('\n Now computing the prediction with a random basis... \n');
Wi = W(knowndata,:);

fprintf('\n QR decomposition of the random basis... \n');;
[Qr,Rr] = qr(Wi);
fprintf('\n Done.\n');

fprintf('\n Solving the least-square problem... \n');
alphar = Rr\(Qr'*knownFy');
fprintf('\n Done.\n');

Fpr = (W*alphar)';

fprintf('\n Computing the quality of the "prediction" (the error should be higher)... \n');
errorpredicrandom = norm(Fpr-Fy)/norm(Fy)
fprintf('\n Done.\n');


end

