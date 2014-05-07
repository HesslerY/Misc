clear

l=8.7;
p=3.7;
h=2.9;

c=3e8;%
f0=3.3e8; %monochromatic pulse frequency

%Reception point coordinates

lambda=c/f0;

X=[(l-lambda)*rand+lambda/2];
Y=[(p-lambda)*rand+lambda/2];
Z=[(h-lambda)*rand+lambda/2];

Nmax=round((l-lambda)*(p-lambda)*(h-lambda)*.74/(4/3*pi*(3e8/f0/4)^3))
while (length(X)<150)

    X_0=(l-lambda)*rand+lambda/2;
    Y_0=(p-lambda)*rand+lambda/2;
    Z_0=(h-lambda)*rand+lambda/2;

    D=sqrt((X_0-X).^2+(Y_0-Y).^2+(Z_0-Z).^2);
    if min(D)>.7*lambda
        X=[X;X_0];
        Y=[Y;Y_0];
        Z=[Z;Z_0];
length(X)
    end
end

save('XYZ150b.mat','X','Y','Z')
