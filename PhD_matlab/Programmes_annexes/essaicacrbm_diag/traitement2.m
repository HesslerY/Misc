        CA_1=[];
        CRBM_1=[];
        CA_5=[];
        CRBM_5=[];
        CA_10=[];
        CRBM_10=[];
        CA_7dB=[];
        CRBM_7dB=[];
        CA_10dB=[];
        CRBM_10dB=[];
        CA_17dB=[];
        CRBM_17dB=[];
        CA_21dB=[];
        CRBM_21dB=[];
        CA1=[];
        CRBM1=[];
        CA5=[];
        CRBM5=[];
        CA10=[];
        CRBM10=[];
        CA7dB=[];
        CRBM7dB=[];
        CA10dB=[];
        CRBM10dB=[];
        CA17dB=[];
        CRBM17dB=[];
        CA21dB=[];
        CRBM21dB=[];

M=30
for q=1:M;
    for j=1:100
        r=[];
        r(1)=randi(149);
        while length(r)<q
            a=randi(149);
            if min(a-r)~=0
                r=[r; a];
            end
        end
        
        CA_1(j)=max(ad1ca(r,60000)/max(ad1ca(:,60000)));
        CRBM_1(j)=max(ad1crbm(r,60000)/max(ad1crbm(:,60000)));
        CA_5(j)=max(ad5ca(r,60000)/max(ad5ca(:,60000)));
        CRBM_5(j)=max(ad5crbm(r,60000)/max(ad5crbm(:,60000)));
        CA_10(j)=max(ad10ca(r,60000)/max(ad10ca(:,60000)));
        CRBM_10(j)=max(ad10crbm(r,60000)/max(ad10crbm(:,60000)));
        CA_7dB(j)=max(ad7dBca(r,60000)/max(ad7dBca(:,60000)));
        CRBM_7dB(j)=max(ad7dBcrbm(r,60000)/max(ad7dBcrbm(:,60000)));
        CA_10dB(j)=max(ad10dBca(r,60000)/max(ad10dBca(:,60000)));
        CRBM_10dB(j)=max(ad10dBcrbm(r,60000)/max(ad10dBcrbm(:,60000)));
        CA_17dB(j)=max(ad17dBca(r,60000)/max(ad17dBca(:,60000)));
        CRBM_17dB(j)=max(ad17dBcrbm(r,60000)/max(ad17dBcrbm(:,60000)));
        CA_21dB(j)=max(ad21dBca(r,60000)/max(ad21dBca(:,60000)));
        CRBM_21dB(j)=max(ad21dBcrbm(r,60000)/max(ad21dBcrbm(:,60000)));
    end
    CA1(q)=mean(CA_1);
    CA5(q)=mean(CA_5);
    CA10(q)=mean(CA_10);
    CA7dB(q)=mean(CA_7dB);
    CA10dB(q)=mean(CA_10dB);
    CA17dB(q)=mean(CA_17dB);
    CA21dB(q)=mean(CA_21dB);
    
    CRBM1(q)=mean(CRBM_1);
    CRBM5(q)=mean(CRBM_5);
    CRBM10(q)=mean(CRBM_10);
    CRBM7dB(q)=mean(CRBM_7dB);
    CRBM10dB(q)=mean(CRBM_10dB);
    CRBM17dB(q)=mean(CRBM_17dB);
    CRBM21dB(q)=mean(CRBM_21dB);
    
    
end

p=1;
figure(1)
subplot(2,2,1)
plot(1:M,filter(ones(p,1)/p,1,CA1),1:M,filter(ones(p,1)/p,1,CRBM1))
ylim([0 1])
subplot(2,2,2)
plot(1:M,filter(ones(p,1)/p,1,CA5),1:M,filter(ones(p,1)/p,1,CRBM5))
ylim([0 1])
subplot(2,2,3)
plot(1:M,filter(ones(p,1)/p,1,CA10),1:M,filter(ones(p,1)/p,1,CRBM10))
ylim([0 1])

figure(2)
subplot(2,2,1)
plot(1:M,filter(ones(p,1)/p,1,CA7dB),1:M,filter(ones(p,1)/p,1,CRBM7dB))
ylim([0 1])
subplot(2,2,2)
plot(1:M,filter(ones(p,1)/p,1,CA10dB),1:M,filter(ones(p,1)/p,1,CRBM10dB))
ylim([0 1])
subplot(2,2,3)
plot(1:M,filter(ones(p,1)/p,1,CA17dB),1:M,filter(ones(p,1)/p,1,CRBM17dB))
ylim([0 1])
subplot(2,2,4)
plot(1:M,filter(ones(p,1)/p,1,CA21dB),1:M,filter(ones(p,1)/p,1,CRBM21dB))
ylim([0 1])