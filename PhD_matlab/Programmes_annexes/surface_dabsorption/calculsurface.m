clear all
i=1;
l=0.61;
h=0.15;
A=[];
f=figure(1);
for n=1:4
    
    voxel([-l/2 -l/2 -n*h/2],[l l n*h],'k',1)
    axis equal
    axis vis3d
    axis off
    xlim([-.8,.8])
    ylim([-.8,.8])
    zlim([-.8,.8])
    for theta=-90:1:90;
        for phi=0:1:359
            view(phi,theta)
            getframe;
            filename=sprintf('fig%d.bmp',i);
            saveas(f,filename,'bmpmono')
            importfile(filename)
            A=[A;n theta phi sum(sum(cdata))/3];
            delete(filename)
            %       close all
            i=i+1;
        end
        
    end
    save('aire.mat','A')
end
