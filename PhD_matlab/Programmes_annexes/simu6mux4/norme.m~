clear all

load FFT8X50.mat

df=f(2)-f(1);

%80-600 MHz
f80=round(80e6/df);
f600=round(600e6/df);

for i=f80:f600
    Emx1(i)=max(abs(Rx1(:,i)));
    Emx2(i)=max(abs(Rx2(:,i)));
    Emx3(i)=max(abs(Rx3(:,i)));
    Emx4(i)=max(abs(Rx4(:,i)));
    Emx5(i)=max(abs(Rx5(:,i)));
    Emx6(i)=max(abs(Rx6(:,i)));
    Emx7(i)=max(abs(Rx7(:,i)));
    Emx8(i)=max(abs(Rx8(:,i)));
    
    Emy1(i)=max(abs(Ry1(:,i)));
    Emy2(i)=max(abs(Ry2(:,i)));
    Emy3(i)=max(abs(Ry3(:,i)));
    Emy4(i)=max(abs(Ry4(:,i)));
    Emy5(i)=max(abs(Ry5(:,i)));
    Emy6(i)=max(abs(Ry6(:,i)));
    Emy7(i)=max(abs(Ry7(:,i)));
    Emy8(i)=max(abs(Ry8(:,i)));

    Emz1(i)=max(abs(Rz1(:,i)));
    Emz2(i)=max(abs(Rz2(:,i)));
    Emz3(i)=max(abs(Rz3(:,i)));
    Emz4(i)=max(abs(Rz4(:,i)));
    Emz5(i)=max(abs(Rz5(:,i)));
    Emz6(i)=max(abs(Rz6(:,i)));
    Emz7(i)=max(abs(Rz7(:,i)));
    Emz8(i)=max(abs(Rz8(:,i)));
end

%600-1200MHz
f1200=round(1200e6/df);

for i=f600+1:f1200
    Emx1(i)=max(abs(Rx1(1:18,i)));
    Emx2(i)=max(abs(Rx2(1:18,i)));
    Emx3(i)=max(abs(Rx3(1:18,i)));
    Emx4(i)=max(abs(Rx4(1:18,i)));
    Emx5(i)=max(abs(Rx5(1:18,i)));
    Emx6(i)=max(abs(Rx6(1:18,i)));
    Emx7(i)=max(abs(Rx7(1:18,i)));
    Emx8(i)=max(abs(Rx8(1:18,i)));
    
    Emy1(i)=max(abs(Ry1(1:18,i)));
    Emy2(i)=max(abs(Ry2(1:18,i)));
    Emy3(i)=max(abs(Ry3(1:18,i)));
    Emy4(i)=max(abs(Ry4(1:18,i)));
    Emy5(i)=max(abs(Ry5(1:18,i)));
    Emy6(i)=max(abs(Ry6(1:18,i)));
    Emy7(i)=max(abs(Ry7(1:18,i)));
    Emy8(i)=max(abs(Ry8(1:18,i)));

    Emz1(i)=max(abs(Rz1(1:18,i)));
    Emz2(i)=max(abs(Rz2(1:18,i)));
    Emz3(i)=max(abs(Rz3(1:18,i)));
    Emz4(i)=max(abs(Rz4(1:18,i)));
    Emz5(i)=max(abs(Rz5(1:18,i)));
    Emz6(i)=max(abs(Rz6(1:18,i)));
    Emz7(i)=max(abs(Rz7(1:18,i)));
    Emz8(i)=max(abs(Rz8(1:18,i)));
end

%1200-2000MHz
f2000=round(2000e6/df);

for i=f1200+1:f2000
    Emx1(i)=max(abs(Rx1(1:12,i)));
    Emx2(i)=max(abs(Rx2(1:12,i)));
    Emx3(i)=max(abs(Rx3(1:12,i)));
    Emx4(i)=max(abs(Rx4(1:12,i)));
    Emx5(i)=max(abs(Rx5(1:12,i)));
    Emx6(i)=max(abs(Rx6(1:12,i)));
    Emx7(i)=max(abs(Rx7(1:12,i)));
    Emx8(i)=max(abs(Rx8(1:12,i)));
    
    Emy1(i)=max(abs(Ry1(1:12,i)));
    Emy2(i)=max(abs(Ry2(1:12,i)));
    Emy3(i)=max(abs(Ry3(1:12,i)));
    Emy4(i)=max(abs(Ry4(1:12,i)));
    Emy5(i)=max(abs(Ry5(1:12,i)));
    Emy6(i)=max(abs(Ry6(1:12,i)));
    Emy7(i)=max(abs(Ry7(1:12,i)));
    Emy8(i)=max(abs(Ry8(1:12,i)));

    Emz1(i)=max(abs(Rz1(1:12,i)));
    Emz2(i)=max(abs(Rz2(1:12,i)));
    Emz3(i)=max(abs(Rz3(1:12,i)));
    Emz4(i)=max(abs(Rz4(1:12,i)));
    Emz5(i)=max(abs(Rz5(1:12,i)));
    Emz6(i)=max(abs(Rz6(1:12,i)));
    Emz7(i)=max(abs(Rz7(1:12,i)));
    Emz8(i)=max(abs(Rz8(1:12,i)));
end

for i=f80:f2000
s8x(i)=std([Emx1(i);Emx2(i);Emx3(i);Emx4(i);Emx5(i);Emx6(i);Emx7(i);Emx8(i)]);
s8y(i)=std([Emy1(i);Emy2(i);Emy3(i);Emy4(i);Emy5(i);Emy6(i);Emy7(i);Emy8(i)]);
s8z(i)=std([Emz1(i);Emz2(i);Emz3(i);Emz4(i);Emz5(i);Emz6(i);Emz7(i);Emz8(i)]);

Emx(i)=mean([Emx1(i);Emx2(i);Emx3(i);Emx4(i);Emx5(i);Emx6(i);Emx7(i);Emx8(i)]);
Emy(i)=mean([Emy1(i);Emy2(i);Emy3(i);Emy4(i);Emy5(i);Emy6(i);Emy7(i);Emy8(i)]);
Emz(i)=mean([Emz1(i);Emz2(i);Emz3(i);Emz4(i);Emz5(i);Emz6(i);Emz7(i);Emz8(i)]);

s24(i)=std([Emx1(i);Emx2(i);Emx3(i);Emx4(i);Emx5(i);Emx6(i);Emx7(i);Emx8(i);Emy1(i);Emy2(i);Emy3(i);Emy4(i);Emy5(i);Emy6(i);Emy7(i);Emy8(i);Emz1(i);Emz2(i);Emz3(i);Emz4(i);Emz5(i);Emz6(i);Emz7(i);Emz8(i)]);
Em24(i)=mean([Emx1(i);Emx2(i);Emx3(i);Emx4(i);Emx5(i);Emx6(i);Emx7(i);Emx8(i);Emy1(i);Emy2(i);Emy3(i);Emy4(i);Emy5(i);Emy6(i);Emy7(i);Emy8(i);Emz1(i);Emz2(i);Emz3(i);Emz4(i);Emz5(i);Emz6(i);Emz7(i);Emz8(i)]);


end


s8xdB=20*log10((s8x+Emx)./Emx);
s8ydB=20*log10((s8y+Emy)./Emy);
s8zdB=20*log10((s8z+Emz)./Emz);

s24dB=20*log10((s24+Em24)./Em24);



gabarit=3*ones(1,f2000);
f400=round(400e6/df);
f100=round(100e6/df);
gabarit(1:f100)=NaN;
gabarit(f100:f400)=linspace(4,3,f400-f100+1);
 
plot(f(1:f2000)/1e6,gabarit,f(1:f2000)/1e6,s8xdB,f(1:f2000)/1e6,s8ydB,f(1:f2000)/1e6,s8zdB)

%entre 80et 600 20fr
for i=1:20
U1(i)=round(80e6*1.106^i/df);
end
%entre 600 et 1200 15fr
for i=1:15
U2(i)=round(600e6*1.0473^i/df);
end
%entre 1200 et 2000 10fr
for i=1:10
U3(i)=round(1200e6*1.0^i/df);
end

U=[U1 U2 U3];


plot(f(1:f2000)/1e6,gabarit,f(U)/1e6,s8xdB(U),f(U)/1e6,s8ydB(U),f(U)/1e6,s8zdB(U))

plot(f(1:f2000)/1e6,gabarit,f(U)/1e6,s24dB(U))


