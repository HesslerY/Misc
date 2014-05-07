for i=1:101
    imagesc((3*lamb),NN(:,1,1),effrand(:,:,i))
    Colorbar;
    caxis([0 100])
    title([num2str(Crit(1,i,1)),' dB'])
    xlabel('lambda')
    ylabel('N')
    getframe;
end


for i=1:101
    imagesc(3*lamb,NN(:,1,1),tempsrand(:,:,i))
    Colorbar;
   caxis([0 200])
    title([num2str(Crit(1,i,1)),' dB'])
    xlabel('lambda')
    ylabel('N')
    getframe; 
end


for i=1:101
    A=effrand(:,100,i);
    B=tempsrand(:,100,i);
    plotyy(1:100,A(:),1:100,B(:))
    title([num2str(Crit(1,i,1)),' dB'])
    grid on
    xlabel('N')
    getframe;
    pause(.2)
end


for i=1:101
    N=10
    A=effrand(N,:,i);
    B=tempsrand(N,:,i);
    plotyy(lamb*3,filter(ones(1,1)/1,1,A(:)),lamb*3,filter(ones(1,1)/1,1,B(:)))
    title([num2str(Crit(1,i,1)),' dB'])
    grid on
    xlabel('lambda')
    getframe;
    pause(.2)
end