load Resultsourcesproches.mat

YY=2:.005:2.3;
for f=2:1:4000
    for i=1:1:61
        A=corrcoef(FFTx(1,:,f),FFTx(i,:,f));
        B=corrcoef(FFTy(1,:,f),FFTy(i,:,f));
        C=corrcoef(FFTz(1,:,f),FFTz(i,:,f));
        Result(i,1)=(YY(i)-YY(1))*freq(f)/3e8;
        Result(i,2)=A(2,1);
        Result(i,3)=B(2,1);
        Result(i,4)=C(2,1);
        
    end

plot(Result(:,1),Result(:,2),Result(:,1),Result(:,3),Result(:,1),Result(:,4))
title([num2str(round(freq(f)/1e6)),' MHz'])
grid on
xlabel('lambda')
ylabel('Corr')
legend('E_x','E_y','E_z')
xlim([0 max(Result(:,1))])
ylim([-0.2 1])
getframe;
filename=sprintf('spatialcorr%d.png',f);
saveas(gcf,filename)
end


YY=2:.005:2.3;
for f=2:10:4000
    for i=1:1:61
        A=corrcoef(FFTx(1,:,f),FFTy(i,:,f));
        B=corrcoef(FFTy(1,:,f),FFTz(i,:,f));
        C=corrcoef(FFTx(1,:,f),FFTz(i,:,f));
        Result(i,1)=(YY(i)-YY(1))*freq(f)/3e8;
        Result(i,2)=A(2,1);
        Result(i,3)=B(2,1);
        Result(i,4)=C(2,1);
        
    end

plot(Result(:,1),Result(:,2),Result(:,1),Result(:,3),Result(:,1),Result(:,4))
title([num2str(round(freq(f)/1e6)),' MHz'])
grid on
xlabel('lambda')
ylabel('Corr')
legend('xy','yz','xz')
xlim([0 max(Result(:,1))])
ylim([-0.2 1])
getframe;
filename=sprintf('componentcorr%d.png',f);
saveas(gcf,filename)
end
