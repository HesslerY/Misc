
clear
tic
%for f0=[4000e6, 3000e6, 2000e6, 1500e6, 1000e6,900e6, 800e6,700e6,600e6,500e6,400e6,300e6,200e6,100e6]
for f0=[200e6]

o1=randi(360);
th1=2*pi*rand;
ph1=pi*rand;

o2=randi(360);
th2=2*pi*rand;
ph2=pi*rand;

o3=randi(360);
th3=2*pi*rand;
ph3=pi*rand;

o4=randi(360);
th4=2*pi*rand;
ph4=pi*rand;

o5=randi(360);
th5=2*pi*rand;
ph5=pi*rand;

o6=randi(360);
th6=2*pi*rand;
ph6=pi*rand;

load freq.mat
frequence=round(f0/freq(2));

q=4;
t=1;


disp('1/6')
u1=0;
for j=1:q:360
    u1=1+u1;
    filex=sprintf('Resultx_%d_%d',1,mod(j+o1,360)+1);
    filey=sprintf('Resulty_%d_%d',1,mod(j+o1,360)+1);
    filez=sprintf('Resultz_%d_%d',1,mod(j+o1,360)+1);

    load(filex)
    load(filey)
    load(filez)

    fftx1(:,u1)=xFFTx(:,frequence)*sin(th1)*cos(ph1+j*2*pi/360);
    ffty1(:,u1)=xFFTy(:,frequence)*sin(th1)*sin(ph1+j*2*pi/360);
    fftz1(:,u1)=xFFTz(:,frequence)*cos(th1);
end
 toc
disp('2/6')
 u2=0;
 for j=1:q:360
    u2=1+u2;
     filex=sprintf('Resultx_%d_%d',2,mod(j+o2,360)+1);
     filey=sprintf('Resulty_%d_%d',2,mod(j+o2,360)+1);
     filez=sprintf('Resultz_%d_%d',2,mod(j+o2,360)+1); 

     load(filex)
     load(filey)
     load(filez)
 
     fftx2(:,u2)=xFFTx(:,frequence)*sin(th2)*cos(ph2+j*2*pi/360);
     ffty2(:,u2)=xFFTy(:,frequence)*sin(th2)*sin(ph2+j*2*pi/360);
     fftz2(:,u2)=xFFTz(:,frequence)*cos(th2);
 end    
 toc
disp('3/6')
u3=0;
for j=1:q:360
    u3=1+u3;
    filex=sprintf('Resultx_%d_%d',3,mod(j+o3,360)+1);
    filey=sprintf('Resulty_%d_%d',3,mod(j+o3,360)+1);
    filez=sprintf('Resultz_%d_%d',3,mod(j+o3,360)+1);

    load(filex)
    load(filey)
    load(filez)

    fftx3(:,u3)=xFFTx(:,frequence)*sin(th3)*cos(ph3+j*2*pi/360);
    ffty3(:,u3)=xFFTy(:,frequence)*sin(th3)*sin(ph3+j*2*pi/360);
    fftz3(:,u3)=xFFTz(:,frequence)*cos(th3);
end    

 toc
disp('4/6')    
 u4=0;
 for j=1:q:360
     u4=1+u4;
     filex=sprintf('Resultx_%d_%d',4,mod(j+o4,360)+1);
     filey=sprintf('Resulty_%d_%d',4,mod(j+o4,360)+1);
     filez=sprintf('Resultz_%d_%d',4,mod(j+o4,360)+1);
 
     load(filex)
     load(filey)
     load(filez)
 
     fftx4(:,u4)=xFFTx(:,frequence)*sin(th4)*cos(ph4+j*2*pi/360);
     ffty4(:,u4)=xFFTy(:,frequence)*sin(th4)*sin(ph4+j*2*pi/360);
     fftz4(:,u4)=xFFTz(:,frequence)*cos(th4);
 end    

 toc
 disp('5/6')
 u5=0;
 for j=1:q:360
     u5=1+u5;
     filex=sprintf('Resultx_%d_%d',5,mod(j+o5,360)+1);
     filey=sprintf('Resulty_%d_%d',5,mod(j+o5,360)+1);
     filez=sprintf('Resultz_%d_%d',5,mod(j+o5,360)+1);
 
     load(filex)
     load(filey)
     load(filez)
 
     fftx5(:,u5)=xFFTx(:,frequence)*sin(th5)*cos(ph5+j*2*pi/360);
     ffty5(:,u5)=xFFTy(:,frequence)*sin(th5)*sin(ph5+j*2*pi/360);
     fftz5(:,u5)=xFFTz(:,frequence)*cos(th5);
 end    
  
 toc   
disp('6/6')
u6=0;
for j=1:q:360
    u6=1+u6;
    filex=sprintf('Resultx_%d_%d',6,mod(j+o6,360)+1);
    filey=sprintf('Resulty_%d_%d',6,mod(j+o6,360)+1);
    filez=sprintf('Resultz_%d_%d',6,mod(j+o6,360)+1);

    load(filex)
    load(filey)
    load(filez)

    fftx6(:,u6)=xFFTx(:,frequence)*sin(th6)*cos(ph6+j*2*pi/360);
    ffty6(:,u6)=xFFTy(:,frequence)*sin(th6)*sin(ph6+j*2*pi/360);
    fftz6(:,u6)=xFFTz(:,frequence)*cos(th6);
end    
 toc
disp('Combination')
fftx=fftx1+fftx2+fftx3+fftx4+fftx5+fftx6;%
ffty=ffty1+ffty2+ffty3+ffty4+ffty5+ffty6;%
fftz=fftz1+fftz2+fftz3+fftz4+fftz5+fftz6;%


for u=1:length(fftx(1,:))
    A=corrcoef(abs(fftx(:,length(fftx(1,:))/2)),abs(fftx(:,u)));
    B=corrcoef(abs(ffty(:,length(ffty(1,:))/2)),abs(ffty(:,u)));
    C=corrcoef(abs(fftz(:,length(fftz(1,:))/2)),abs(fftz(:,u)));
    Rx(u,t)=A(2,1);
    Ry(u,t)=B(2,1);
    Rz(u,t)=C(2,1);
end
toc
phi=0:2*pi/length(fftx(1,:)):2*pi-2*pi/length(fftx(1,:));
plot(phi*360/2/pi,mean(Rx,2),phi*360/2/pi,mean(Ry,2),phi*360/2/pi,mean(Rz,2))
xlabel('°')
ylabel('Corr')
title([num2str(round(freq(frequence)/1e6)),' MHz'])
legend('E_x','E_y','E_z')
grid on
xlim([0 360])
ylim([-.2 1.1])
% 
% save('Result6morceaux.mat','phi1','Rx','Ry','Rz','freq','frequence')
%filename=sprintf('5palesgraph_%d.fig',round(freq(frequence)/1e6))
%filename2=sprintf('5palesresult_%d.mat',round(freq(frequence)/1e6))
%saveas(gcf,filename)
%save(filename2,'Rx','Ry','Rz','phi')
end
