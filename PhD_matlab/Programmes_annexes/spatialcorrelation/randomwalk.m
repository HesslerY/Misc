l=8.7
p=3.7
h=2.9

X=[1]
Y=[1]
Z=[1]

N=500
while(length(X)<N)
    a=length(X)
    pr=(8.7+2.9+3.7)*rand;
    if pr<2.9
        X=[X;X(a)];
            Y=[Y;Y(a)+0.015];
            Z=[Z;Z(a)];
    else
    if pr<(2.9+2)
         X=[X;X(a)];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)+0.015];
    else
     
         X=[X;X(a)+0.015];
            Y=[Y;Y(a)];
            Z=[Z;Z(a)]; 
     end
    end
end

max([X Y Z])