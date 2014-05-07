clear

E=[];
c=3e8;
D=3;
for j=[0,4,5,6]
    s=pi*(.1)^2*j
    R=.998;
    for i=1:3000
        Nbray(i)=4*i^2+2;
        absorbant(i)=(4*i^2+2)*s;
        surf(i)=4*pi*i^2*D^2;
        dist(i)=i*D;
        ordre(i)=R^i;
    end
    P=1-absorbant./surf;
    E=[E;cumprod(P).*ordre];
    
end
plot((1:3000)*D/c,E)
xlim([0 30e-6])
