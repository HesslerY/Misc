load result4_6MU.mat

Sx=Signalfinalx;
Sz=Signalfinalz;

load result3_6MU.mat

Sx(1:113,:)=Signalfinalx;
Sz(1:113,:)=Signalfinalz;

load result2_6MU.mat

Sx(1:75,:)=Signalfinalx;
Sz(1:75,:)=Signalfinalz;

load result1_6MU.mat

Sx(1:38,:)=Signalfinalx;
Sz(1:38,:)=Signalfinalz;


for i=1:150
    Sxenvp(i,:)=filter(ones(100,1)/100,1,Sx(i,:).^2);
    Szenvp(i,:)=filter(ones(100,1)/100,1,Sz(i,:).^2);
end