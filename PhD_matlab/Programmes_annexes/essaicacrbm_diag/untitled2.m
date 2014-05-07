for i=1:150
plot(tt,CONVant(i,:),'.',tt,CONVomni(i,:),'o')
getframe;
pause
end

p=100;
for i=1:150
    omnf=filter(ones(p,1)/p,1,abs(CONVomni10dB(i,:)));
    omn(i)=omnf(60000);
       antf=filter(ones(p,1)/p,1,abs(CONVant10dB(i,:)));
    ant(i)=antf(60000);
end

mean(omn)
mean(ant)


plot(1:150,omn,1:150,ant)

% 
% plot(1:150,A,1:150,B)