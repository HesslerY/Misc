%  x=1-(0:.0001:1)
%  y=log(x)/log(1-exp(-0.025))
% 
% plot(x,y*360)



for i=1:10000
A=10*log10(exprnd(1,360,1));

Q=find(A>8);
R=find(A>7);
S=find(A>6);
T=find(A>5);
U=find(A>4);
V=find(A>3);
W=find(A>2);
X=find(A>1);
Y=find(A>0);
m(i)=max(A);
QQ(i)=length(Q);
RR(i)=length(R);
SS(i)=length(S);
TT(i)=length(T);
UU(i)=length(U);
VV(i)=length(V);
WW(i)=length(W);
XX(i)=length(X);
YY(i)=length(Y);
end

Pmax=mean(m)
psup8dB=mean(QQ)/360
psup7dB=mean(RR)/360
psup6dB=mean(SS)/360
psup5dB=mean(TT)/360
psup4dB=mean(UU)/360
psup3dB=mean(VV)/360
psup2dB=mean(WW)/360
psup1dB=mean(XX)/360
psup0dB=mean(YY)/360



largeursup4dB=mean(UU)/30
largeursup5dB=mean(TT)/22
largeursup6dB=mean(SS)/15


f=-5:.1:10;
g=360*log(0.95)./log(1-exp(-1*10.^(-f./10)));

plot(f,g)
%largeur à 3dB???



c=10.^(f/10);
h=(1-exp(-c))./(log(0.91)./log(1-exp(-c)));
plot(f,h)

B=exprnd(1,10000,1)
Y=sort(B)
