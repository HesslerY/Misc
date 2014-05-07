D=2.8141;

Rcopo=0.9981;
Rxpo=0.9983;
c=299792458;
%traitement 1-4 absorbant

importfiles2p('mesures_aire/1surblocxpo.s2p')
surblocxpo1=data;
t=data(:,1);
surblocxpo1=complex(surblocxpo1(:,4),surblocxpo1(:,5));
Psurblocxpo1=cumsum(abs(surblocxpo1).^2)/sum(abs(surblocxpo1).^2);
U=find(Psurblocxpo1>.63);
tau_1=t(U(1));


S1=(pi*D^2)*(1-exp(-D/(2*c*tau_1))/Rxpo);


importfiles2p('mesures_aire/2surblocxpo.s2p')
surblocxpo2=data;
surblocxpo2=complex(surblocxpo2(:,4),surblocxpo2(:,5));
Psurblocxpo2=cumsum(abs(surblocxpo2).^2)/sum(abs(surblocxpo2).^2);
U=find(Psurblocxpo2>.63);
tau_2=t(U(1));


S2=(pi*D^2)*(1-exp(-D/(2*c*tau_2))/Rxpo);

importfiles2p('mesures_aire/3surblocxpo.s2p')
surblocxpo3=data;
surblocxpo3=complex(surblocxpo3(:,4),surblocxpo3(:,5));
Psurblocxpo3=cumsum(abs(surblocxpo3).^2)/sum(abs(surblocxpo3).^2);
U=find(Psurblocxpo3>.63);
tau_3=t(U(1));


S3=(pi*D^2)*(1-exp(-D/(2*c*tau_3))/Rxpo);

importfiles2p('mesures_aire/4surblocxpo.s2p')
surblocxpo4=data;
surblocxpo4=complex(surblocxpo4(:,4),surblocxpo4(:,5));
Psurblocxpo4=cumsum(abs(surblocxpo4).^2)/sum(abs(surblocxpo4).^2);
U=find(Psurblocxpo4>.63);
tau_4=t(U(1));


S4=(pi*D^2)*(1-exp(-D/(2*c*tau_4))/Rxpo);



%pyramide sur bloc

importfiles2p('mesures_aire/pyramidebloc.s2p')
pysurbloc=data;
pysurbloc=complex(pysurbloc(:,4),pysurbloc(:,5));
Ppysurbloc=cumsum(abs(pysurbloc).^2)/sum(abs(pysurbloc).^2);
U=find(Ppysurbloc>.63);
tau_pybloc=t(U(1));


S_pybloc=(pi*D^2)*(1-exp(-D/(2*c*tau_pybloc))/Rxpo);


%pyramide sur sol

importfiles2p('mesures_aire/pyramidesol.s2p')
pysursol=data;
pysursol=complex(pysursol(:,4),pysursol(:,5));
Ppysursol=cumsum(abs(pysursol).^2)/sum(abs(pysursol).^2);
U=find(Ppysursol>.63);
tau_pysol=t(U(1));


S_pysol=(pi*D^2)*(1-exp(-D/(2*c*tau_pysol))/Rxpo);


