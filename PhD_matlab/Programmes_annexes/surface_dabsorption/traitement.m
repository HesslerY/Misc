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




%traitement 1 absorbant

importfiles2p('mesures_aire/1deboutxpo.s2p')
deboutxpo=data;
deboutxpo=complex(deboutxpo(:,4),deboutxpo(:,5));
Pdeboutxpo=cumsum(abs(deboutxpo).^2)/sum(abs(deboutxpo).^2);
U=find(Pdeboutxpo>.63);
tau_deboutxpo=t(U(1));

Sdeboutxpo=(pi*D^2)*(1-exp(-D/(2*c*tau_deboutxpo))/Rxpo);



importfiles2p('mesures_aire/1arretexpo.s2p')
arretexpo=data;
arretexpo=complex(arretexpo(:,4),arretexpo(:,5));
Parretexpo=cumsum(abs(arretexpo).^2)/sum(abs(arretexpo).^2);
U=find(Parretexpo>.63);
tau_arretexpo=t(U(1));

Sarretexpo=(pi*D^2)*(1-exp(-D/(2*c*tau_arretexpo))/Rxpo);


importfiles2p('mesures_aire/1coinxpo.s2p')
coinxpo=data;
coinxpo=complex(coinxpo(:,4),coinxpo(:,5));
Pcoinxpo=cumsum(abs(coinxpo).^2)/sum(abs(coinxpo).^2);
U=find(Pcoinxpo>.63);
tau_coinxpo=t(U(1));

Scoinxpo=(pi*D^2)*(1-exp(-D/(2*c*tau_coinxpo))/Rxpo);


importfiles2p('mesures_aire/1couchexpo.s2p')
couchexpo=data;
couchexpo=complex(couchexpo(:,4),couchexpo(:,5));
Pcouchexpo=cumsum(abs(couchexpo).^2)/sum(abs(couchexpo).^2);
U=find(Pcouchexpo>.63);
tau_couchexpo=t(U(1));

Scouchexpo=(pi*D^2)*(1-exp(-D/(2*c*tau_couchexpo))/Rxpo);

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




%pyramide sur bloc 3x4

importfiles2p('mesures_aire/py34.s2p')
pysurbloc34=data;
pysurbloc34=complex(pysurbloc34(:,4),pysurbloc34(:,5));
Ppysurbloc34=cumsum(abs(pysurbloc34).^2)/sum(abs(pysurbloc34).^2);
U=find(Ppysurbloc34>.63);
tau_pybloc34=t(U(1));


S_pybloc34=(pi*D^2)*(1-exp(-D/(2*c*tau_pybloc34))/Rxpo);