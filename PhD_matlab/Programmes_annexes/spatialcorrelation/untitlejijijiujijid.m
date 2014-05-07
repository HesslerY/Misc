clear
load ResultHill500elements_x.mat
load ResultHill500elements_y.mat
load ResultHill500elements_z.mat
load ResultHill500elements_111.mat


%load ResultHill4x100elementsb.mat
xFFTx=xFFTx(:,1:10000).^2;
xFFTy=xFFTy(:,1:10000).^2;
xFFTz=xFFTz(:,1:10000).^2;

yFFTx=yFFTx(:,1:10000).^2;
yFFTy=yFFTy(:,1:10000).^2;
yFFTz=yFFTz(:,1:10000).^2;

zFFTx=zFFTx(:,1:10000).^2;
zFFTy=zFFTy(:,1:10000).^2;
zFFTz=zFFTz(:,1:10000).^2;

s111FFTx=s111FFTx(:,1:10000).^2;
s111FFTy=s111FFTy(:,1:10000).^2;
s111FFTz=s111FFTz(:,1:10000).^2;


xFFTx(1:60,:)=[];
xFFTy(1:60,:)=[];
xFFTz(1:60,:)=[];

yFFTx(1:60,:)=[];
yFFTy(1:60,:)=[];
yFFTz(1:60,:)=[];

zFFTx(1:60,:)=[];
zFFTy(1:60,:)=[];
zFFTz(1:60,:)=[];

s111FFTx(1:60,:)=[];
s111FFTy(1:60,:)=[];
s111FFTz(1:60,:)=[];

figure(1)
k=0;
% for f=1405:1:5244 
%for f=2:1400;%[524,1049,2097,4194,6291,8389] 
for f=316:1574;%[524,1049,2097,4194,6291,8389] 
k=k+1
    for p=1:61
        A=corrcoef(xFFTx(1:61:5978,f),xFFTx(p:61:5978,f)); %l
        B=corrcoef(xFFTy(1:61:5978,f),xFFTy(p:61:5978,f)); %t
        C=corrcoef(xFFTz(1:61:5978,f),xFFTz(p:61:5978,f)); %t

        D=corrcoef(yFFTx(1:61:5978,f),yFFTx(p:61:5978,f)); %t
        E=corrcoef(yFFTy(1:61:5978,f),yFFTy(p:61:5978,f)); %l
        F=corrcoef(yFFTz(1:61:5978,f),yFFTz(p:61:5978,f)); %t

        G=corrcoef(zFFTx(1:61:5978,f),zFFTx(p:61:5978,f)); %t
        H=corrcoef(zFFTy(1:61:5978,f),zFFTy(p:61:5978,f)); %t
        I=corrcoef(zFFTz(1:61:5978,f),zFFTz(p:61:5978,f)); %l

        J=corrcoef(s111FFTx(1:61:5978,f),s111FFTx(p:61:5978,f)); %n
        K=corrcoef(s111FFTy(1:61:5978,f),s111FFTy(p:61:5978,f)); %n
        L=corrcoef(s111FFTz(1:61:5978,f),s111FFTz(p:61:5978,f)); %n


        Resultl(1,p)=A(2,1);
        Resultl(2,p)=E(2,1);
        Resultl(3,p)=I(2,1);

        Resultt(1,p)=B(2,1);
        Resultt(2,p)=C(2,1);
        Resultt(3,p)=D(2,1);
        Resultt(4,p)=F(2,1);
        Resultt(5,p)=G(2,1);
        Resultt(6,p)=H(2,1);

        Resultn(1,p)=J(2,1);
        Resultn(2,p)=K(2,1);
        Resultn(3,p)=L(2,1);


    end
%     for p=1:61
%         Result(p)=0.005*(p-1)*2*pi*freq(f)/3e8;
%     end

   p=1:61;
         Result=0.005*(p-1).*2*pi*freq(f)/3e8;
     %endfor i=1:length(Result(:,1))

     for i=1:61
         d=Result(1,i);
         
         Resulth(i,2)=((3/2*(sin(d)/d-1/d^2*(sin(d)/d-cos(d)))))^2;
         Resulth(i,1)=(3/d^2*(sin(d)/d-cos(d)))^2;
         Resulth(i,3)=(sin(d)/d)^2;
     end

    %subplot(2,2,k)
    plot(Result,mean(Resultl),'-b',Result,Resulth(:,1),'--b',Result,mean(Resultt),'-r',Result,Resulth(:,2),'--r',Result,mean(Resultn),'-k',Result,Resulth(:,3),'--k')
    grid on
    xlim([0 max(Result)])
    ylim([-0.2 1.05])
    %h=legend('$\rho_{ll}$ sim.','$\rho_{ll}$ th.','$\rho_{tt}$ sim.','$\rho_{tt}$ th.','$\rho^2$ sim.','$\rho^2$ th.')
    %set(h,'Interpreter','Latex')
    title([num2str(round(freq(f)/1e6)),' MHz'],'Interpreter','Latex','FontSize',14)
    xlabel('$kr$','Interpreter','Latex','FontSize',12)
    ylabel('Corr.','Interpreter','Latex','FontSize',12)

    getframe;
       filename=sprintf('figd_%d.png',f)
       saveas(gcf,filename)

end
