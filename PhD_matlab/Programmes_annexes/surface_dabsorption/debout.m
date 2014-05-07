%debout
     
     
clear all
close all


h=.15
l=.60
w=.60


fig=figure(1)

% voxel([-l/2 -h/2 0],[l,h,0.001],[1,1,1],1)
% voxel([-l/2 -h/2 w],[l,h,0.001],[0,0,0],1)
% voxel([-l/2 -h/2 0],[l,0.001,w],[0,0,0],1)
% voxel([-l/2 h/2 0],[l,0.001,w],[0,0,0],1)
% voxel([l/2 -h/2 0],[0.001,h,w],[0,0,0],1)
% voxel([-l/2 -h/2 0],[0.001,h,w],[0,0,0],1)
voxel([-l/2 -w/2 0],[l,w,0.001],[0,0,0],1)
voxel([-l/2 -w/2 h],[l,w,0.001],[0,0,0],1)
voxel([-l/2 -w/2 0],[l,0.001,h],[0,0,0],1)
voxel([-l/2 w/2 0],[l,0.001,h],[0,0,0],1)
voxel([l/2 -w/2 0],[0.001,w,h],[1,1,1],1)
voxel([-l/2 -w/2 0],[0.001,w,h],[0,0,0],1)
axis equal
axis vis3d
axis off

    xlim([-.8,.8])
    ylim([-.8,.8])
    zlim([0,.8])
Ap=[];
i=0;
    for theta=-90:90;
        for phi=0:1:359
            view(phi,theta)
            pause(0.4)
            getframe;
            filename=sprintf('figd%d.bmp',i);
            saveas(fig,filename,'bmpmono')
            importfile(filename)
            Ap=[Ap;theta phi sum(sum(cdata))/3];
            delete(filename)
            %       close all
            i=i+1;
        end
        
    end
     save('debout.mat','Ap')
     
     
 % Define these variables appropriately:
mail = 'manuamador@gmail.com'; %Your GMail email address
password = 'kashmir3534'; %Your GMail password

% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('manuamador@gmail.com','debout','debout','debout.mat')
