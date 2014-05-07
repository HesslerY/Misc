%load Result.mat
F=[];

for f=1:1:1100
    Ax=[];
    Ay=[];
    Az=[];
    for i=1:50;

        Ax=[Ax;FFTx(i,:,f)];
        Ay=[Ay;FFTy(i,:,f)];
        Az=[Az;FFTz(i,:,f)];
    end

    decox=[Ax(:,1)];
    decoy=[Ay(:,1)];
    decoz=[Az(:,1)];
    for i=1:210
        bx=[];
        for j=1:length(decox(1,:))

            CCx=corrcoef(Ax(:,i),decox(:,j));
            bx(j)=abs(CCx(2,1));
        end
        if max(bx(j))<0.15
            decox=[decox Ax(:,i)];
        end
        
        
              by=[];
        for j=1:length(decoy(1,:))

            CCy=corrcoef(Ay(:,i),decoy(:,j));
            by(j)=abs(CCy(2,1));
        end
        if max(by(j))<0.15
            decoy=[decoy Ay(:,i)];
        end
        
        
              bz=[];
        for j=1:length(decoz(1,:))

            CCz=corrcoef(Az(:,i),decoz(:,j));
            bz(j)=abs(CCz(2,1));
        end
        if max(bz(j))<0.15
            decoz=[decoz Az(:,i)];
        end
    end
    
    F=[F;freq(f)/1e6 length(decox(1,:)) length(decoy(1,:)) length(decoz(1,:))]
end
for o=1:1100
    Nmax(o)=round(1*.74/(4/3*pi*(3e8/freq(o)/4)^3))
end

plot(F(:,1),filter(ones(1,50)/50,1,F(:,2)),F(:,1),filter(ones(1,50)/50,1,F(:,3)),F(:,1),filter(ones(1,50)/50,1,F(:,4)),F(:,1),Nmax)
xlim([0 1000])
grid on



