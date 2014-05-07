load Resulttemporapido.mat
Sx=filter(ones(50,1)/50,1,abs(Signalx));
Sy=filter(ones(50,1)/50,1,abs(Signaly));
Sz=filter(ones(50,1)/50,1,abs(Signalz));

% for i=2:10:20000
% hist(mean(Sx(1:i,:)'),50)
% title([num2str(t(i)/1e-9),' ns'])
% getframe;
% end
% for i=1:50:20000
%     %subplot(3,1,1)
%     hist(Sx(i,:),20)
%     title([num2str(t(i)/1e-9),' ns'])
%     getframe;
% end
%
%
% for i=1:100
%     %subplot(3,1,1)
%     plot(t,Sx(:,i))
%     %title([num2str(t(i)/1e-9),' ns'])
%     pause(1)
%     getframe;
% end
p=5000
for i=1:100:20000-p
    for j=1:p
        A= corrcoef(Sx(i,:),Sx(i+j,:));
        B=  corrcoef(Sy(i,:),Sy(i+j,:));
        C=corrcoef(Sz(i,:),Sz(i+j,:));
        R(i+j,1)=t(i+j)/1e-9;
        R(i+j,2)=A(2,1);
        R(i+j,3)=B(2,1);
        R(i+j,4)=C(2,1);
    end


end

plot(R(:,1),R(:,2),R(:,1),R(:,3),R(:,1),R(:,4))

for i=1:1:20000
Rayl(i,1)=t(i)/1e-9;
Rayl(i,2)=raylfit(Sx(i,:));
Rayl(i,3)=raylfit(Sy(i,:));
Rayl(i,4)=raylfit(Sz(i,:));
Rayl(i,5)=std(Sx(i,:))/mean(Sx(i,:));
Rayl(i,6)=std(Sy(i,:))/mean(Sy(i,:));
Rayl(i,7)=std(Sz(i,:))/mean(Sz(i,:));
Rayl(i,8)=max(Sx(i,:));
Rayl(i,9)=max(Sy(i,:));
Rayl(i,10)=max(Sz(i,:));
end

plot(Rayl(:,1),Rayl(:,2),Rayl(:,1),Rayl(:,3),Rayl(:,1),Rayl(:,4))%,Rayl(:,1),Rayl(:,5),Rayl(:,1),Rayl(:,6),Rayl(:,1),Rayl(:,7))
grid on

plot(Rayl(:,1),Rayl(:,5),Rayl(:,1),Rayl(:,6),Rayl(:,1),Rayl(:,7))
grid on

plot(Rayl(:,1),Rayl(:,8),Rayl(:,1),Rayl(:,9),Rayl(:,1),Rayl(:,10))
grid on



for i=1:1:19999
    drayl(i,1)=t(i)/1e-9;
    drayl(i,2)=(Rayl(i+1,2)-Rayl(i,2))/(t(2)-t(1));
    drayl(i,3)=(Rayl(i+1,3)-Rayl(i,3))/(t(2)-t(1));
    drayl(i,4)=(Rayl(i+1,4)-Rayl(i,4))/(t(2)-t(1));
end
plot(drayl(:,1),drayl(:,2),drayl(:,1),drayl(:,3),drayl(:,1),drayl(:,4))
grid on



for i=500:1:20000
Rayl(i,1)=t(i)/1e-9;
A=wblfit(Sx(i,:));
B=wblfit(Sy(i,:));
C=wblfit(Sz(i,:));

Rayl(i,2)=A(2);
Rayl(i,3)=B(2);
Rayl(i,4)=C(2);
end

plot(Rayl(:,1),Rayl(:,2),Rayl(:,1),Rayl(:,3),Rayl(:,1),Rayl(:,4))













