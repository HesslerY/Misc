YY=2:.005:2.3;
for f=1749
    for i=1:1:61
        A=corrcoef(FFTx(1,:,f),FFTx(i,:,f));
        B=corrcoef(FFTy(1,:,f),FFTy(i,:,f));
        C=corrcoef(FFTz(1,:,f),FFTz(i,:,f));
        Result1000(i,1)=(YY(i)-YY(1))*freq(f)/3e8
        Result1000(i,2)=A(2,1);
        Result1000(i,3)=B(2,1);
        Result1000(i,4)=C(2,1);
        
    end
end

YY=2:.005:2.3;
for f=875
    for i=1:1:61
        A=corrcoef(FFTx(1,:,f),FFTx(i,:,f));
        B=corrcoef(FFTy(1,:,f),FFTy(i,:,f));
        C=corrcoef(FFTz(1,:,f),FFTz(i,:,f));
        Result500(i,1)=(YY(i)-YY(1))*freq(f)/3e8
        Result500(i,2)=A(2,1);
        Result500(i,3)=B(2,1);
        Result500(i,4)=C(2,1);
        
    end
end

YY=2:.005:2.3;
for f=2624
    for i=1:1:61
        A=corrcoef(FFTx(1,:,f),FFTx(i,:,f));
        B=corrcoef(FFTy(1,:,f),FFTy(i,:,f));
        C=corrcoef(FFTz(1,:,f),FFTz(i,:,f));
        Result1500(i,1)=(YY(i)-YY(1))*freq(f)/3e8
        Result1500(i,2)=A(2,1);
        Result1500(i,3)=B(2,1);
        Result1500(i,4)=C(2,1);
        
    end
end

YY=2:.005:2.3;
for f=3498
    for i=1:1:61
        A=corrcoef(FFTx(1,:,f),FFTx(i,:,f));
        B=corrcoef(FFTy(1,:,f),FFTy(i,:,f));
        C=corrcoef(FFTz(1,:,f),FFTz(i,:,f));
        Result2000(i,1)=(YY(i)-YY(1))*freq(f)/3e8
        Result2000(i,2)=A(2,1);
        Result2000(i,3)=B(2,1);
        Result2000(i,4)=C(2,1);
        
    end
end









subplot(2,2,1)
plot(Result500(:,1),Result500(:,2),Result500(:,1),Result500(:,3),Result500(:,1),Result500(:,4))
title('500 MHz')
grid on
xlabel('lambda')
ylabel('Corr')
legend('E_x','E_y','E_z')
xlim([0 .5])
subplot(2,2,2)
plot(Result1000(:,1),Result1000(:,2),Result1000(:,1),Result1000(:,3),Result1000(:,1),Result1000(:,4))
title('1 GHz')
grid on
xlabel('lambda')
ylabel('Corr')
legend('E_x','E_y','E_z')
xlim([0 1])

subplot(2,2,3)
plot(Result1500(:,1),Result1500(:,2),Result1500(:,1),Result1500(:,3),Result1500(:,1),Result1500(:,4))
title('1500 MHz')
grid on
xlabel('lambda')
ylabel('Corr')
legend('E_x','E_y','E_z')
xlim([0 1.5])

subplot(2,2,4)
plot(Result2000(:,1),Result2000(:,2),Result2000(:,1),Result2000(:,3),Result2000(:,1),Result2000(:,4))
title('2 GHz')
grid on
xlabel('lambda')
ylabel('Corr')
legend('E_x','E_y','E_z')
xlim([0 2])


