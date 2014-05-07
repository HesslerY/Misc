clear
Result=[];
RICE=[];
ParamV=[];
s=0.08;
dim=100;
Nexp=1000

vvv=-31:.2:2

for v=10.^(vvv/10)

    x = s .* randn(dim,1) + v;
    y = s .* randn(dim,1);
    r = sqrt(x.^2 + y.^2);
    RICE=[RICE r];
    ParamV=[ParamV;v];
end
tic
alpha=45;
uuu=-30:1:1;

for u=10.^(uuu/10)
    %u=[0.0005 .001 .002 .003 .004 .005 .006 .007 .009 .01 .015 .02 .03 .05 .08 .1 .15 .2 .3 .4 .5 .6]
    u
    RR1=[];
    for k=1:1:50

        a=(100*k/50);
        if mod(a,20)==0
            disp(a)
        end

        v=abs(complex(0.0339268*randn(Nexp,dim)+u*sind(45),0.0339268*randn(Nexp,dim)+u*cosd(45)));
        sigma=0.0085*randn(Nexp,dim)+0.08;
        S=abs(complex(sigma.*randn(Nexp,dim)+(v*cosd(45)),sigma.*randn(Nexp,dim)+(v*sind(45))));
        N=length(S(1,:)); %posiiton de brasseur ou nombre d'�chantillons

      
        D=[];%matrice des valeurs de d
        M=0;
        for f=1:1:length(S(:,1)) %Boucle sur les fr�quences
            R=abs(S(f,:)); % on prend le module des valeurs de la fr�quence consid�r�e
            for j=1:1:length(RICE(1,:))
                [Verdict(f,j), PROB(f,j)]=kstest2(R,RICE(:,j));
            end


        end
    end

for i=1:1:length(PROB(:,1))
    V(i)=(PROB(i,:))*ParamV;
    %if V(i)~=0 
        V(i)=V(i)/sum(PROB(i,:));
    %end
end
Result=[Result;u/(sqrt(2)*s) mean(V)/(sqrt(2)*s)]
end
toc
%save ResultN100-10000.mat 'Result'