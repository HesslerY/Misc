%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                               %%%%%%
%%%%%       CHANNEL IMPULSE RESPONSE, PULSED SIGNAL RESPONSE        %%%%%%
%%%%%                    & FREQUENCY RESPONSE                       %%%%%%
%%%%%              by E. Amador (manuamador@gmail.com)              %%%%%%
%%%%%                         IETR/DGA                              %%%%%%
%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear

freq=[];
l=8.7;
p=3.7;
h=2.9;

c=3e8;%
Lt=1e-6; %Time-window length in seconds
nbre_elements=1
dmax=c*Lt %maximal distance

tau=.3e-6; %length of the pulse in seconds
f0=1e9; %monochromatic pulse frequency

N=20000 %number of points for the chosen time-window (Lt)
t=0:Lt/(N-1):Lt; %time
x=0:1/((N-1)/Lt):tau;
s=sin(2*pi*f0*x); %pulsed signal


%Attenuation coefficient for each direction
Rx=0.998;
Ry=0.998;
Rz=0.998;


%Reception point coordinates

load XYZ150b.mat


for v=[6,3]
    disp(v)
    for w=1:360
        disp([num2str(w),'/360'])
        tic
        filename = sprintf('x_%d_%d.mat',v,w);
        load(filename)
        for u=1:1:length(X)


            disp([num2str(100*u/(length(X)))])
            X_1=X(u);
            Y_1=Y(u);
            Z_1=Z(u);

            Sx=zeros(1,N);
            Sy=zeros(1,N);
            Sz=zeros(1,N);
            for j=1:1:length(POS) %loop over the image-currents... to be vctorized... one day may be...

                DX=X_1-POS(j,1);
                DY=Y_1-POS(j,2);
                DZ=Z_1-POS(j,3);

                dist=sqrt(DX^2+DY^2+DZ^2);
                zl=round((N-1)*dist/c/Lt);

                if zl<N
                    alpha=POS(j,9);
                    beta=POS(j,10);
                    ca=cos(alpha);
                    sa=sin(alpha);
                    cb=cos(beta);
                    sb=sin(beta);

                    %Ralpha/beta, coordinates transformation matrix:
                    %                           Ralpbeta=[  (-sin(beta))^2+(1-(-sin(beta))^2)*cos(alpha)   -sin(beta)*cos(beta)*(1-cos(alpha)) cos(beta)*sin(alpha);
                    %                                     -sin(beta)*cos(beta)*(1-cos(alpha))      (cos(beta))^2+(1-(cos(beta))^2)*cos(alpha) sin(beta)*sin(alpha);
                    %                                     -cos(beta)*sin(alpha)                   -sin(beta)*sin(alpha)                        cos(alpha)];

                    %rectangular coordinates calculation in the local system attached to the considered current (developped expressions, a matrix product is slower)
                    distx=((-sb)^2+(1-(-sb)^2)*ca)*DX+(-sb*cb*(1-ca))*DY+(cb*sa)*DZ;
                    disty=(-sb*cb*(1-ca))*DX+((cb)^2+(1-cb^2)*ca)*DY+(sb*sa)*DZ;
                    distz=(-cb*sa)*DX+(-sb*sa)*DY+(ca)*DZ;

                    distxy=sqrt(distx^2+disty^2);

                    R=Rx^POS(j,4)*Ry^POS(j,5)*Rz^POS(j,6); %Attenuation
                    E=POS(j,8)/dist; % Free-space attenuation
                    costheta=distz/dist;
                    sintheta=distxy/dist;
                    cosphi=distx/distxy;
                    sinphi=disty/distxy;
                    Antth=-sintheta; %dipole radiation pattern

                    %Reverse transformation matrix
                    %                         Ralpbetainv=[  ((-sb)^2+(1-(-sb)^2)*ca) (-sb*cb*(1-ca)) (-cb*sa);
                    %                                     (-sb*cb*(1-ca)) ((cb)^2+(1-(cb)^2)*ca)  (-sb*sa);
                    %                                     (cb*sa) (sb*sa) ca];


                    %Projection in the usual rectangular coordinates

                    Vx=(((-sb)^2+(1-(-sb)^2)*ca)*(Antth*costheta*cosphi)+(-sb*cb*(1-ca))*(Antth*costheta*sinphi)+(-cb*sa)*(-sintheta*Antth));
                    Vy=((-sb*cb*(1-ca))*(Antth*costheta*cosphi)+((cb)^2+(1-(cb)^2)*ca)*(Antth*costheta*sinphi)+(-sb*sa)*(-sintheta*Antth));
                    Vz=((cb*sa)*(Antth*costheta*cosphi)+(sb*sa)*(Antth*costheta*sinphi)+ca*(-sintheta*Antth));


                    %Three axis channel impulse response construction
                    Sx(zl+1)=Sx(zl+1)+R*E*Vx;
                    Sy(zl+1)=Sy(zl+1)+R*E*Vy;
                    Sz(zl+1)=Sz(zl+1)+R*E*Vz;

                end
            end

            %Frequency response
            %disp('FFT')


            Fs = N/Lt; %sampling frequency
            T = 1/Fs;                     % sample length
            L = N;                  %Number of points
            tt = (0:L-1)*T;                % time...

            NFFT = 2^nextpow2(N); % Next power of 2 from length of y
            Yx = fft((Sx),NFFT)/N;
            Yy = fft((Sy),NFFT)/N;
            Yz = fft((Sz),NFFT)/N;
            f = Fs/2*linspace(0,1,NFFT/2);

            % Three axis frequency response
            xFFTx(u,:)=abs(Yx(1:7000));
            xFFTy(u,:)=abs(Yy(1:7000));
            xFFTz(u,:)=abs(Yz(1:7000));
            freq=f(1:7000);

        end
        filename2=sprintf('Resultx_%d_%d.mat',v,w);
        save(filename2,'X','Y','Z','xFFTx','xFFTy','xFFTz','freq','f0')
        toc
        clear xFFTx xFFTy xFFTz
    end
end
