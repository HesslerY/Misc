clear
load ResultHill3x100elements.mat
FFTx=FFTx(:,1:120000).^2;
FFTy=FFTy(:,1:120000).^2;
FFTz=FFTz(:,1:120000).^2;
q=1
j=0;
l=400;


    for f=1:q:length(FFTx(1,:))
        j=j+1;
        A=corrcoef(FFTx(1:4:l,f),FFTx(2:4:l,f));
        B=corrcoef(FFTy(1:4:l,f),FFTy(2:4:l,f));
        C=corrcoef(FFTz(1:4:l,f),FFTz(2:4:l,f));
        
        Result(j,1)=0.015*2*pi*freq(f)/3e8;
        Result_l(j,1)=A(2,1);
        Result_t(j,1)=B(2,1);
        Result_t(j,2)=C(2,1);

        
        D=corrcoef(FFTx(1:4:l,f),FFTx(3:4:l,f));
        E=corrcoef(FFTy(1:4:l,f),FFTy(3:4:l,f));
        F=corrcoef(FFTz(1:4:l,f),FFTz(3:4:l,f));
        

        Result_l(j,2)=E(2,1);
        Result_t(j,3)=D(2,1);
        Result_t(j,4)=F(2,1);
        
        G=corrcoef(FFTx(1:4:l,f),FFTx(4:4:l,f));
        H=corrcoef(FFTy(1:4:l,f),FFTy(4:4:l,f));
        I=corrcoef(FFTz(1:4:l,f),FFTz(4:4:l,f));
        

        Result_l(j,3)=I(2,1);
        Result_t(j,5)=G(2,1);
        Result_t(j,6)=H(2,1);

    end
% for i=1:length(Result(:,1))
% d=Result(i,1);
% Result(i,5)=(3/2*(sin(d/d-1/d^2*(sin(d)/d-cos(d)))));
% Result(i,6)=(1/3*2*Result(i,5)+3/d^2*(sin(d)/d-cos(d)));
% end
% 

p=300
plot(Result(:,1),filter(ones(1,p)/p,1,mean(Result_l,2)),Result(:,1),filter(ones(1,p)/p,1,mean(Result_t,2)))
legend('\rho_{ll}','\rho_{tt}')
grid on
xlabel('kr')
ylabel('Corr')
 xlim([0.1 10])
 ylim([-.2 1])
 