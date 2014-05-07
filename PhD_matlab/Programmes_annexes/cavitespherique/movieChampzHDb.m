% clear
%  load CARTOVHDezgwen200HD.mat
%  load CARTOHHDezgwen200HD.mat
% load CARTOTHDez.mat

%load CARTOHHD.mat
%figure1 = figure('PaperSize',[20 28])

%colormap('hot')
cmax=[];
set(gcf,'nextplot','replacechildren');
x=-.15:0.005:.15;
y=-.15:.005:.15;
z=-.15:.005:.15;
for j=1:1:12000%[50,100,500,1000,2000,3000,5000,7000,10000]
    figure(1)
    
     %subplot(2,3,1)
    %imagesc(x,y,10*log10(CARTOfxH(:,:,j).^2'))
    %caxis([-100 nanmax(nanmax(10*log10(CARTOfxH(:,:,j).^2')))]);
    % title(['Plan H, Ex, ',num2str(ceil(f(j)/1e6)),' MHz'])
    % colorbar('location','Eastoutside')
    % xlabel('x')
    % ylabel('y')
     %axis tight
    % axis image
    %   subplot(2,3,2)
    % imagesc(x,y,10*log10(CARTOfyH(:,:,j).^2'))
     %caxis([-100 nanmax(nanmax(10*log10(CARTOfyH(:,:,j).^2')))]);
    % title(['Plan H, Ey, ',num2str(ceil(f(j)/1e6)),' MHz'])
    % colorbar('location','Eastoutside')
   %  xlabel('x')
   %  ylabel('y')
   %  axis tight
    % axis image
    %  subplot(2,3,3)
      subplot(2,2,4)
    imagesc(x,y,10*log10(CARTOfzH(:,:,j).^2'))
%caxis([-100 max(max((10*log10(CARTOfzH(:,:,j).^2'))))]);
    title(['Plan H, Ez, ',num2str(ceil(f(j)/1e6)),' MHz'])
    colorbar('location','Eastoutside')
    xlabel('x')
    ylabel('y')
    axis tight
    axis image
     
        %subplot(2,3,4)
        subplot(2,2,1)
    imagesc(x,y,10*log10(CARTOfxV(:,:,j).^2'))
%caxis([-100 -30]);
caxis([-100 max(max((10*log10(CARTOfxV(:,:,j).^2'))))]); 
title(['Plan V, Ex, ',num2str(ceil(f(j)/1e6)),' MHz'])
    colorbar('location','Eastoutside')
    xlabel('x')
    ylabel('z')
    axis tight
    axis image
     % subplot(2,3,5)
     subplot(2,2,2)
    imagesc(x,y,10*log10(CARTOfyV(:,:,j).^2'))
  caxis([-100 max(max((10*log10(CARTOfyV(:,:,j).^2'))))]); 
    % caxis([-100 0]);
    title(['Plan V, Ey, ',num2str(ceil(f(j)/1e6)),' MHz'])
    colorbar('location','Eastoutside')
    xlabel('x')
    ylabel('z')
    axis tight
    axis image
      %subplot(2,3,6)
      subplot(2,2,3)
    imagesc(x,y,10*log10(CARTOfzV(:,:,j).^2'))
    %caxis([-100 max(max(10*log10(CARTOfzV(:,:,j).^2')))]);
    title(['Plan V, Ez, ',num2str(ceil(f(j)/1e6)),' MHz'])
    colorbar('location','Eastoutside')
    xlabel('x')
    ylabel('z')
    axis tight
    axis image
    
    
getframe;
          filename = sprintf('sphere_gwenn200dB2x2HD_%d.png',j);
% % %          %filename2 = sprintf('sscolorbarcartofreqH_%d.pdf',j);
       saveas(gcf,filename)%,'epsc')
%       
% %      %saveas(gca,filename2)



end
