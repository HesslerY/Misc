%reproduction resultat de Hill


clear

global Lt c Rx Ry Rz N POS

freq=[];
l=8.7;
p=3.7;
h=2.9;

c=3e8;%
Lt=1e-6; %Time-window length in seconds
nbre_elements=1
dmax=c*Lt %maximal distance

tau=1e-6; %length of the pulse in seconds
f0=1e9; %monochromatic pulse frequency

N=round(20*Lt*f0) %number of points for the chosen time-window (Lt)
t=0:Lt/(N-1):Lt; %time
x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %pulsed signal


%Attenuation coefficient for each direction
Rx=0.96;
Ry=0.96;
Rz=0.96;


%Reception point coordinates

lambda=c/f0;

X=[(l-lambda)*rand+lambda/2];
Y=[(p-lambda)*rand+lambda/2];
Z=[(h-lambda)*rand+lambda/2];
Nmax=round((l-lambda)*(p-lambda)*(h-lambda)*.74/(4/3*pi*(3e8/f0/4)^3))
while (length(X)<100)
    
    X_0=(l-lambda)*rand+lambda/2;
    Y_0=(p-lambda)*rand+lambda/2;
    Z_0=(h-lambda)*rand+lambda/2;
    
    D=sqrt((X_0-X).^2+(Y_0-Y).^2+(Z_0-Z).^2);
    if min(D)>lambda/2
        X=[X;X_0];
        Y=[Y;Y_0];
        Z=[Z;Z_0];
    end
end


    
    load('1elem_1000ns7.mat')
    for u=1:length(X)
        tic
        disp([num2str(u),'/100'])
        X_1=X(u);
        Y_1=Y(u);
        Z_1=Z(u);
        
       % disp('CIR')
        [Sx,Sy,Sz]=CIR(X_1,Y_1,Z_1);
     sigX=conv(Sx,s);
     sigY=conv(Sy,s);
     sigZ=conv(Sz,s);
     
    Signalx(:,u)=sigX(1:N);
    Signaly(:,u)=sigY(1:N);
    Signalz(:,u)=sigZ(1:N);
    toc
    end
    
toc
save('Resulttemporapido.mat','Signalx','Signaly','Signalz','X','Y','Z','t','f0')
