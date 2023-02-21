function EOF(datafile,doplot,threshold)
% datafile = name of the file containing data or 'gendata' to use the synthetic data generator
% doplot = true if plots are requested
% threshold = threshold for filtering eigenvalues, in terms of percentage of variability to be recovered

close all;

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
 addnoise=0;
 mydata=gendata(nx,ny,nt,neof,0);
 lon=1:nx;
 lat=1:ny;
 time=1:nt;
end

nx=length(lon);
ny=length(lat);
nt=length(time);

%% 1/ EOF on the whole data minus last year.
F=mydata(1:nt-12,:);

% Remove undefined values
F(F==-32768)=0;

% Animated plot (can be very long, use with care...)
if(doplot)
  figure; 
  fprintf('ENTER to go on\n');
  pause;
  fprintf('Displaying...\n');
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
      stopplot=input('0 to stop, ENTER to go on.\n');
      if(stopplot==0)
        break;
      end
    end
  end
  fprintf('ENTER to go on\n');
  pause;
end
end
