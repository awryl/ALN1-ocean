function mydata = gendata(nx,ny,nt,neof,addnoise)
% Generates synthetic data.
% Input:
%  nx: size of the domain on the x-axis
%  ny: size of the domain on the y-axis
%  nt: number of time steps
%  neof: number of non-zeros eigenvalues in the spectrum
%  addnoise: true to add some random noise
% Ouput:
%  mydata: matrix of size nt times (nx*ny).
%
% The output data is such it yields neof EOFs. Each EOF is
% a 1x1 square (i.e. 4 points on the grid), i.e. is non-zero
% on that square and 0 outside,  and is associated to a
% cosine principal components. The EOFs are disjoint spatially
% (squares dont touch one another) and the associated PCs
% have different amplitudes and periods. The spectrum
% of the covariance matrix should exhibit clearly separated
% eigenvalues.
%
% Adding noise should modify the quality of the EOF basis
% and the quality of the prediction. The noise is quite
% limited in amplitude so that extracted EOFs should 
% remain the same as without noise.
%
% NB: nx and ny have to be larger than 4*neof so that there
% is enough space on the domain to generate disjoint EOFs.

% Check conditions
if(nx<4*neof | ny<4*neof)
 fprintf('nx>=4*neof and ny>=4*neof must hold.\n');
 mydata=[];
 return
end

% Amplitudes: so that eigenvalues are fairly separated
alpha=1+0.5*(1:neof)';

% Periods: so that we have enough information/oscillations 
omega=40*(1:neof)';

% Principal components: each PC_i is alpha(i)cos(omega(i)*t)
PC=cos((1/nt:1/nt:1)'*omega')*diag(alpha);

% Empirical orthogonal functions: each EOF_i is 1 inside 2x2 square
% (defined by a random corner) near the center of the grid, 0 outside

 % Lower-left corner
 llcx = floor(nx/4)+(0:floor(nx/2/neof):floor(nx/2/neof)*neof);
 llcy = floor(ny/4)+(0:floor(ny/2/neof):floor(ny/2/neof)*neof);

 % Upper-right corner
 urcx = llcx+1;
 urcy = llcy+1;

 % EOF
 EOF=zeros(neof,nx*ny);
 for i=1:neof
  tmp=zeros(nx,ny);
  tmp(llcx(i):urcx(i),llcy(i):urcy(i))=1;
  EOF(i,:)=reshape(tmp,nx*ny,1);
 end

% Data: sum_i PC_i EOF_i
mydata=PC*EOF;

% Optionally, add some random noise.
if(addnoise)
 mydata=mydata+rand(nt,nx*ny)/min(min(mydata))/10;
end

end
