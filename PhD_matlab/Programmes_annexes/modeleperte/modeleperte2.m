clear

E=[];
c=3e8;
D=3;
for j=0:6
    s=pi*(.096)^2*j %surface efficace (disque) * nombre d'objets
    R=.9981; %perte chambre vide (parois, antennes...)
    for i=1:3000
        absorbant(i)=(4*i^2+2)*s; %surface d'absorbants pour l'ordre i
        surf(i)=4*pi*i^2*D^2;   %surface de la sphère à l'ordre i
        ordre(i)=R^i; %atténuation des parois pour un ordre i
    end
    P=1-(absorbant./surf); %probabilité  pour un rayon de passer à travers l'ordre i
    A=(cumprod(P).*ordre); %probabilité combinée pour un rayon de passer à travers les i ordres * atténuation de l'ordre i
    E=[E;A/max(A)]; %normalisation du champ
    
end
plot((1:3000)*D/c/1e-9,E)
xlim([0 30000])
grid on
xlabel('Time in ns','Interpreter','Latex')
ylabel('$\vert E_{norm} \vert$','Interpreter','Latex')
h=legend('0 abs','1 abs','2 abs','3 abs','4 abs','5 abs','6 abs');
set(h,'Interpreter','Latex')