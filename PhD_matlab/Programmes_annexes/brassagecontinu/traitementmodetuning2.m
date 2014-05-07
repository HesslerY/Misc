u=0
for i=[21,51,71,101]
    u=1+u
    subplot(2,2,u)
    imagesc(NN(:,1,1),(3*lamb*100),effrand(:,:,i)')
    %Colorbar;
    caxis([0 100])
    title(['Niveau: ',num2str(Crit(1,i,1)),' dB'])
    xlabel('$N$')
   ylabel('$\lambda$ [cm]')
    getframe;
end

u=0
for i=[21,51,71,101]
    u=1+u
    subplot(2,2,u)
    imagesc(NN(:,1,1),(3*lamb*100),tempsrand(:,:,i)'/360*90)
    Colorbar;
   %caxis([0 100])
    title(['Niveau: ',num2str(Crit(1,i,1)),' dB'])
    xlabel('N')
    ylabel('$\lambda$ [cm]')
    getframe; 
end


u=0
for i=[21,51,71,101]
    u=1+u
    subplot(2,2,u)
    A=effrand(:,1000,i);
    B=tempsrand(:,1000,i)/360*90;
    plotyy(1:100,A(:),1:100,B(:))
    title(['Niveau: ',num2str(Crit(1,i,1)),' dB'])
    grid on
    xlabel('$N$')
    getframe;
    pause(.2)
end

p=10
u=0
for i=[21,51,71,101]
    u=1+u
    subplot(2,2,u)
    N=100
    A=effrand(N,:,i);
    B=tempsrand(N,:,i)/360*90;
    plotyy(lamb*3*100,filter(ones(1,p)/p,1,A(:)),lamb*3*100,filter(ones(1,p)/p,1,B(:)))
    title(['Niveau: ',num2str(Crit(1,i,1)),' dB'])
    grid on
    xlabel('$\lambda$ [cm]')
    getframe;
    pause(.2)
end