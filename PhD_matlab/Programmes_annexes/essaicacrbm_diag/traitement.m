clear all
load essaiCRBMD10dB.mat
CRBMa10dB=CONVant10dB;
CRBMo10dB=CONVomni10dB;
CRBMa10dB(81,:)=[];
CRBMo10dB(81,:)=[];
clear CONVant10dB CONVomni10dB


load essaiCRBMD7dB.mat
CRBMa7dB=CONVant7dB;
CRBMo7dB=CONVomni7dB;
CRBMa7dB(81,:)=[];
CRBMo7dB(81,:)=[];
clear CONVant7dB CONVomni7dB

load essaiCRBMD17dB.mat
CRBMa17dB=CONVant17dB;
CRBMo17dB=CONVomni17dB;
CRBMa17dB(81,:)=[];
CRBMo17dB(81,:)=[];
clear CONVant17dB CONVomni17dB

load essaiCRBMD21dB.mat
CRBMa21dB=CONVant21dB;
CRBMo21dB=CONVomni21dB;
CRBMa21dB(81,:)=[];
CRBMo21dB(81,:)=[];
clear CONVant21dB CONVomni21dB

load essaiCRBMD1n5.mat
CRBMa1=CONVant1;
CRBMo1=CONVomni1;
CRBMa1(81,:)=[];
CRBMo1(81,:)=[];
clear CONVant1 CONVomni1

load essaiCRBMD5n5.mat
CRBMa5=CONVant5;
CRBMo5=CONVomni5;
CRBMa5(81,:)=[];
CRBMo5(81,:)=[];
clear CONVant5 CONVomni5

load essaiCRBMD10n5.mat
CRBMa10=CONVant10;
CRBMo10=CONVomni10;
CRBMa10(81,:)=[];
CRBMo10(81,:)=[];
clear CONVant10 CONVomni10

%%%%%%%%%%%%%%%%%%%

load essaiCAD10dB.mat
CAa10dB=CONVant;
CAo10dB=CONVomni;
CAa10dB(81,:)=[];
CAo10dB(81,:)=[];
clear CONVant CONVomni


load essaiCAD7dB.mat
CAa7dB=CONVant;
CAo7dB=CONVomni;
CAa7dB(81,:)=[];
CAo7dB(81,:)=[];
clear CONVant CONVomni



load essaiCAD17dB.mat
CAa17dB=CONVant;
CAo17dB=CONVomni;
CAa17dB(81,:)=[];
CAo17dB(81,:)=[];
clear CONVant CONVomni



load essaiCAD21dB.mat
CAa21dB=CONVant;
CAo21dB=CONVomni;
CAa21dB(81,:)=[];
CAo21dB(81,:)=[];
clear CONVant CONVomni


load essaiCAD1n5.mat
CAa1=CONVant;
CAo1=CONVomni;
CAa1(81,:)=[];
CAo1(81,:)=[];
clear CONVant CONVomni

load essaiCAD5n5.mat
CAa5=CONVant;
CAo5=CONVomni;
CAa5(81,:)=[];
CAo5(81,:)=[];
clear CONVant CONVomni

load essaiCAD10n5.mat
CAa10=CONVant;
CAo10=CONVomni;
CAa10(81,:)=[];
CAo10(81,:)=[];
clear CONVant CONVomni alp bet

p=100;
for i=1:149
    omnf=filter(ones(p,1)/p,1,abs(CAo1(i,:)));
    cao1(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa1(i,:)));
    caa1(i)=antf(60000);
    ao1ca(i,:)=omnf;
    ad1ca(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CAo5(i,:)));
    cao5(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa5(i,:)));
    caa5(i)=antf(60000);
    
    ao5ca(i,:)=omnf;
    ad5ca(i,:)=antf;

    
    omnf=filter(ones(p,1)/p,1,abs(CAo10(i,:)));
    cao10(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa10(i,:)));
    caa10(i)=antf(60000);
    
    ao10ca(i,:)=omnf;
    ad10ca(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CAo10dB(i,:)));
    cao10dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa10dB(i,:)));
    caa10dB(i)=antf(60000);
    
    ao10dBca(i,:)=omnf;
    ad10dBca(i,:)=antf;
     
    
    omnf=filter(ones(p,1)/p,1,abs(CAo7dB(i,:)));
    cao7dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa7dB(i,:)));
    caa7dB(i)=antf(60000);
    
    ao7dBca(i,:)=omnf;
    ad7dBca(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CAo17dB(i,:)));
    cao17dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa17dB(i,:)));
    caa17dB(i)=antf(60000);
    
    ao17dBca(i,:)=omnf;
    ad17dBca(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CAo21dB(i,:)));
    cao21dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CAa21dB(i,:)));
    caa21dB(i)=antf(60000);
    
    ao21dBca(i,:)=omnf;
    ad21dBca(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CRBMo1(i,:)));
    cro1(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa1(i,:)));
    cr1(i)=antf(60000);
    
    ao1crbm(i,:)=omnf;
    ad1crbm(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CRBMo5(i,:)));
    cro5(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa5(i,:)));
    cr5(i)=antf(60000);
    
    ao5crbm(i,:)=omnf;
    ad5crbm(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CRBMo10(i,:)));
    cro10(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa10(i,:)));
    cr10(i)=antf(60000);
    
    ao10crbm(i,:)=omnf;
    ad10crbm(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CRBMo10dB(i,:)));
    cro10dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa10dB(i,:)));
    cr10dB(i)=antf(60000);
    
    ao10dBcrbm(i,:)=omnf;
    ad10dBcrbm(i,:)=antf;
   
    omnf=filter(ones(p,1)/p,1,abs(CRBMo7dB(i,:)));
    cro7dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa7dB(i,:)));
    cr7dB(i)=antf(60000);
    
    ao7dBcrbm(i,:)=omnf;
    ad7dBcrbm(i,:)=antf;
    
    omnf=filter(ones(p,1)/p,1,abs(CRBMo17dB(i,:)));
    cro17dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa17dB(i,:)));
    cr17dB(i)=antf(60000);
    
    ao17dBcrbm(i,:)=omnf;
    ad17dBcrbm(i,:)=antf;
        
    omnf=filter(ones(p,1)/p,1,abs(CRBMo21dB(i,:)));
    cro21dB(i)=omnf(60000);
    antf=filter(ones(p,1)/p,1,abs(CRBMa21dB(i,:)));
    cr21dB(i)=antf(60000);
    
    ao21dBcrbm(i,:)=omnf;
    ad21dBcrbm(i,:)=antf;
    
end


figure(1)
plot(tt,10*log10(max(ad1ca.^2)./mean(ad1ca.^2)))%,tt,

figure(2)
plot(tt,10*log10(max(ad5ca.^2)./mean(ad5ca.^2)))%,tt,

figure(3)
plot(tt,10*log10(max(ad10ca.^2)./mean(ad10ca.^2)))%,tt,

figure(4)
plot(tt,10*log10(max(ad10dBca.^2)./mean(ad10dBca.^2)))%,tt,


figure(1)
plot(tt,10*log10(max(ad1crbm.^2)./mean(ad1crbm.^2)))%,tt,

figure(2)
plot(tt,10*log10(max(ad5crbm.^2)./mean(ad5crbm.^2)))%,tt,

figure(3)
plot(tt,10*log10(max(ad10crbm.^2)./mean(ad10crbm.^2)))%,tt,

figure(4)
plot(tt,10*log10(max(ad10dBcrbm.^2)./mean(ad10dBcrbm.^2)))%,tt,

figure(1)
plot(tt,10*log10(max(ao1crbm.^2)./mean(ao1crbm.^2)))%,tt,

figure(2)
plot(tt,10*log10(max(ao5crbm.^2)./mean(ao5crbm.^2)))%,tt,

figure(3)
plot(tt,10*log10(max(ao10crbm.^2)./mean(ao10crbm.^2)))%,tt,

figure(4)
plot(tt,10*log10(max(ao10dBcrbm.^2)./mean(ao10dBcrbm.^2)))%,tt,


plot(tt,mean(ao1crbm.^2)/mean(mean(ao1crbm.^2)),tt,mean(ad1crbm.^2)/mean(mean(ad1crbm.^2)),tt,mean(ad5crbm.^2)/mean(mean(ad5crbm.^2)),tt,mean(ad10crbm.^2)/mean(mean(ad10crbm.^2)),tt,mean(ad10dBcrbm.^2)/mean(mean(ad10dBcrbm.^2)))

plot(tt,mean(ao1ca.^2)/mean(mean(ao1ca.^2)),tt,mean(ad1ca.^2)/mean(mean(ad1ca.^2)),tt,mean(ad5ca.^2)/mean(mean(ad5ca.^2)),tt,mean(ad10ca.^2)/mean(mean(ad10ca.^2)),tt,mean(ad10dBca.^2)/mean(mean(ad10dBca.^2)))

subplot(4,1,1)
plot(1:150,caa1,1:150,caa5,1:150,caa10),
subplot(4,1,2)
plot(1:150,cr1,1:150,cr5,1:150,cr10)
subplot(4,1,3)
plot(1:150,caa7dB,1:150,caa10dB,1:150,caa17dB,1:150,caa21dB)
subplot(4,1,4)
plot(1:150,cr7dB,1:150,cr10dB,1:150,cr17dB,1:150,cr21dB)

figure(4)
plot(tt,var(ao1crbm./mean(mean(ao1crbm))),tt,var(ad1crbm./mean(mean(ad1crbm))),tt,var(ad5crbm./mean(mean(ad5crbm))),tt,var(ad10crbm./mean(mean(ad10crbm))),tt,var(ad10dBcrbm./mean(mean(ad10dBcrbm))))%,tt,

figure(5)
plot(tt,var(ao1ca./mean(mean(ao1ca))),tt,var(ad1ca/mean(mean(ad1ca))),tt,var(ad5ca./mean(mean(ad5ca))),tt,var(ad10ca./mean(mean(ad10ca))),tt,var(ad10dBca./mean(mean(ad10dBca))))%,tt,


plot(1:150,ao1crbm(:,60000),1:150,ad1crbm(:,60000))
plot(1:150,ao5crbm(:,60000),1:150,ad5crbm(:,60000))
plot(1:150,ao10crbm(:,60000),1:150,ad10crbm(:,60000))
plot(1:150,ao10dBcrbm(:,60000),1:150,ad10dBcrbm(:,60000))



plot(1:150,ao1ca(:,60000),1:150,ad1ca(:,60000))
plot(1:150,ao5ca(:,60000),1:150,ad5ca(:,60000))
plot(1:150,ao10ca(:,60000),1:150,ad10ca(:,60000))
plot(1:150,ao10dBca(:,60000),1:150,ad10dBca(:,60000))

q=3

figure(1)
plot(1:q:150,ad1ca(1:q:150,60000)/max(ad1ca(:,60000)),1:q:150,ad1crbm(1:q:150,60000)/max(ad1crbm(:,60000)))
figure(2)
plot(1:q:150,ad5ca(1:q:150,60000)/max(ad5ca(:,60000)),1:q:150,ad5crbm(1:q:150,60000)/max(ad5crbm(:,60000)))
figure(3)
plot(1:q:150,ad10ca(1:q:150,60000)/max(ad10ca(:,60000)),1:q:150,ad10crbm(1:q:150,60000)/max(ad10crbm(:,60000)))
figure(4)
plot(1:q:150,ad10dBca(1:q:150,60000)/max(ad10dBca(:,60000)),1:q:150,ad10dBcrbm(1:q:150,60000)/max(ad10dBcrbm(:,60000)))


for q=1:75
CA1(q)=max(ad1ca(1:q:150,60000)/max(ad1ca(:,60000)))
CRBM1(q)=max(ad1crbm(1:q:150,60000)/max(ad1crbm(:,60000)))
CA5(q)=max(ad5ca(1:q:150,60000)/max(ad5ca(:,60000)))
CRBM5(q)=max(ad5crbm(1:q:150,60000)/max(ad5crbm(:,60000)))
CA10(q)=max(ad10ca(1:q:150,60000)/max(ad10ca(:,60000)))
CRBM10(q)=max(ad10crbm(1:q:150,60000)/max(ad10crbm(:,60000)))
CA7dB(q)=max(ad7dBca(1:q:150,60000)/max(ad7dBca(:,60000)))
CRBM7dB(q)=max(ad7dBcrbm(1:q:150,60000)/max(ad7dBcrbm(:,60000)))
CA10dB(q)=max(ad10dBca(1:q:150,60000)/max(ad10dBca(:,60000)))
CRBM10dB(q)=max(ad10dBcrbm(1:q:150,60000)/max(ad10dBcrbm(:,60000)))
CA17dB(q)=max(ad17dBca(1:q:150,60000)/max(ad17dBca(:,60000)))
CRBM17dB(q)=max(ad17dBcrbm(1:q:150,60000)/max(ad17dBcrbm(:,60000)))
CA21dB(q)=max(ad21dBca(1:q:150,60000)/max(ad21dBca(:,60000)))
CRBM21dB(q)=max(ad21dBcrbm(1:q:150,60000)/max(ad21dBcrbm(:,60000)))
end

p=5;
figure(1)
subplot(2,2,1)
plot(1:75,filter(ones(p,1)/p,1,CA1),1:75,filter(ones(p,1)/p,1,CRBM1))
ylim([0 1])
subplot(2,2,2)
plot(1:75,filter(ones(p,1)/p,1,CA5),1:75,filter(ones(p,1)/p,1,CRBM5))
ylim([0 1])
subplot(2,2,3)
plot(1:75,filter(ones(p,1)/p,1,CA10),1:75,filter(ones(p,1)/p,1,CRBM10))
ylim([0 1])

figure(2)
subplot(2,2,1)
plot(1:75,filter(ones(p,1)/p,1,CA7dB),1:75,filter(ones(p,1)/p,1,CRBM7dB))
ylim([0 1])
subplot(2,2,2)
plot(1:75,filter(ones(p,1)/p,1,CA10dB),1:75,filter(ones(p,1)/p,1,CRBM10dB))
ylim([0 1])
subplot(2,2,3)
plot(1:75,filter(ones(p,1)/p,1,CA17dB),1:75,filter(ones(p,1)/p,1,CRBM17dB))
ylim([0 1])
subplot(2,2,4)
plot(1:75,filter(ones(p,1)/p,1,CA21dB),1:75,filter(ones(p,1)/p,1,CRBM21dB))
ylim([0 1])