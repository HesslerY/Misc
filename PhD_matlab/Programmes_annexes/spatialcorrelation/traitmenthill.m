clear
load ResultHillb.mat
q=10
j=0;
for f=1:q:80000
    j=j+1;

        A=corrcoef(FFTx(1,f:f+q),FFTx(2,f:f+q));
        B=corrcoef(FFTy(1,f:f+q),FFTy(2,f:f+q));
        C=corrcoef(FFTz(1,f:f+q),FFTz(2,f:f+q));
        Result(j,1)=0.015*2*pi*freq(f+q/2)/3e8;
        Result(j,2)=A(2,1);
        Result(j,3)=B(2,1);
        Result(j,4)=C(2,1);

end
for i=1:length(Result(:,1))
d=Result(i,1);
Result(i,5)=(3/2*(sin(d/d-1/d^2*(sin(d)/d-cos(d)))));
Result(i,6)=(1/3*2*Result(i,5)+3/d^2*(sin(d)/d-cos(d)));
end





p=50
plot(Result(:,1),filter(ones(1,p)/p,1,Result(:,2)),Result(:,1),filter(ones(1,p)/p,1,Result(:,3)),Result(:,1),filter(ones(1,p)/p,1,Result(:,4)))%,Result(:,1),Result(:,5),Result(:,1),Result(:,6))


%title([num2str(round(freq(f)/1e6)),' MHz'])
legend('E_x','E_y','E_z')
grid on
xlabel('kr')
ylabel('Corr')
 xlim([0.1 6])
 ylim([-.2 1])
% % 


clear
load ResultHillb.mat
q=100
j=0;
for f=1:q:80000
    j=j+1;

        A=corrcoef(sqrt(FFTx(1,f:f+q).^2+FFTy(1,f:f+q).^2+FFTz(1,f:f+q).^2),sqrt(FFTx(2,f:f+q).^2+FFTy(2,f:f+q).^2+FFTz(2,f:f+q).^2));
   
        Result(j,1)=0.015*2*pi*freq(f+q/2)/3e8;
        Result(j,2)=A(2,1);


end
p=10
plot(Result(:,1),filter(ones(1,p)/p,1,Result(:,2)))%,Result(:,1),filter(ones(1,p)/p,1,Result(:,3)),Result(:,1),filter(ones(1,p)/p,1,Result(:,4)))
%title([num2str(round(freq(f)/1e6)),' MHz'])
legend('E')%,'E_y','E_z')
grid on
xlabel('kr')
ylabel('Corr')
xlim([0.1 6])
ylim([-.2 1])
