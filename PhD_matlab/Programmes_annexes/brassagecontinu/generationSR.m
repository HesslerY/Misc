SR=[]
for i=1:50
X=[];
disp(i)

lambda=1;
N=900000;
theta=2*pi/N:2*pi/N:2*pi;

Npos=360;
thetap=2*pi/Npos:2*pi/Npos:2*pi;

Np=i*2; %Nombre approximatif de pics pos indépendantes du brasseur ?????
while length(X)<N
r=randi([round(N/Np/10) round(N/Np)]);
A=exprnd(lambda,r,1);
A=sort(A);
B=exprnd(lambda,r,1);
B=flipud(sort(B,1));

X=[X;A;B];
end
X=X(1:N);

SR(:,i)=X;
end

save('SR50Np100-2.mat','SR','theta')


q=round(length(X)/Npos)+1;
Xp=X(1:q:length(X));
%Y=circshift(X,randi(length(X)))

%Z=max([X Y]')';


plot(thetap-2*pi/Npos,10*log10(Xp),theta,10*log10(X))





dfittool(X)

dfittool(Xp)