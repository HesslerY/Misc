%empilement
load aire.mat

U=find(A(:,1)==1);
V=find(A(:,1)==2);
W=find(A(:,1)==3);
X=find(A(:,1)==4);

Amin=min(A(U,4));
ratio=0.61*.15/Amin;

A1=mean(A(U,4))*ratio;
A2=mean(A(V,4))*ratio;
A3=mean(A(W,4))*ratio;
A4=mean(A(X,4))*ratio;



%pyramides

load airepyramide.mat

ratiop=0.61/A(1,3);

%pyramides en l'air

Py=mean(A(:,3))*ratiop

%pyramides posées
load plan.mat
Y=find(Ap(:,1)<0);
sol=sum(Ap(Y,3));
Pyp=(sum(A(:,3))-sol)/length(A(:,1))*ratiop


%bloc sur le sol
load aire.mat
A1sol=(sum(A(U,4))-sol)/length(A(U,1))*ratio
A2sol=(sum(A(V,4))-sol)/length(A(V,1))*ratio
A3sol=(sum(A(W,4))-sol)/length(A(W,1))*ratio
A4sol=(sum(A(X,4))-sol)/length(A(X,1))*ratio

A1sol/A1
%2/3
A2sol/A2
%3/4
A3sol/A3
%4/5
A4sol/A4
%5/6




n=4;


l=0.6;
p=0.15;
h=0.6;

R=1/2*sqrt((l^2+p^2)+h^2);

%Omega=1-4*asin(l/R*p/R/sqrt((4*(h/2/R)^2+(l/R)^2)*(4*(h/2/R)^2+(p/R)^2)))/4/pi
Om=1-4*asin(sin(atan(l/h))*sin(atan(p/h)))/4/pi



