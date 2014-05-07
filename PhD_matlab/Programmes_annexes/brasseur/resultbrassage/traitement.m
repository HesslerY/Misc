load optim200.mat
corrm=mean(optim(:,21:23)')';

optim(:,24)=corrm;

optim=sortrows(optim,24);

o1=optim(:,3)*pi/180;
th1=optim(:,4);
phi1=optim(:,5);

o2=optim(:,6)*pi/180;
th2=optim(:,7);
phi2=optim(:,8);

o3=optim(:,9)*pi/180;
th3=optim(:,10);
phi3=optim(:,11);

o4=optim(:,12)*pi/180;
th4=optim(:,13);
phi4=optim(:,14);

o5=optim(:,15)*pi/180;
th5=optim(:,16);
phi5=optim(:,17);

o6=optim(:,18)*pi/180;
th6=optim(:,19);
phi6=optim(:,20);

corrx=optim(:,21);
corry=optim(:,22);
corrz=optim(:,23);

o=[o1 o2 o3 o4 o5 o6]';
th=[th1 th2 th3 th4 th5 th6]';
ph=[phi1 phi2 phi3 phi4 phi5 phi6]';


dth6=abs(th6-th5);
dth5=abs(th5-th4);
dth4=abs(th4-th3);
dth3=abs(th3-th2);
dth2=abs(th2-th1);

dth=mean([dth2 dth3 dth4 dth5 dth6]')';

do6=abs(o6-o5);
do5=abs(o5-o4);
do4=abs(o4-o3);
do3=abs(o3-o2);
do2=abs(o2-o1);

do=mean([do2 do3 do4 do5 do6]')';
corrm=optim(:,24);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot(th1,corrm,'.',th2,corrm,'.',th3,corrm,'.',th4,corrm,'.',th5,corrm,'.',th6,corrm,'.')

plot(dth,corrm)


for i=300:-1:1
hold on
plot(o(:,i),corrm(i),'.','MarkerSize',32)
xlim([0 2*pi])
ylim([0 1])
getframe;
end



