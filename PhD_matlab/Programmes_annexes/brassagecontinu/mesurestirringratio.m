C=[];
M=360*25;

%--Pilotage du DSO

    %--Création d'une interface objet
%dso= tcpip('10.1.45.252');
%set(analyseur, 'OutputBufferSize', 32000');
%set(analyseur, 'InputBufferSize', 32000');
%set(analyseur,'Timeout',10);
%fopen(analyseur);
   dso = visa('ni','TCPIP::10.1.45.252::INSTR');
    set(dso, 'OutputBufferSize', 200000000);
    set(dso, 'InputBufferSize', 20000000);

    %--Ouverture de la connexion réseau DSO

    fopen(dso);

    %--Caractéristiques

    points = str2num(query(dso,'HORizontal:RESOlution?'));
    temps_par_division =str2num( query(dso,'HORizontal:SCAle?'));
    duree=10*temps_par_division;
    t=(10*temps_par_division)/points:(10*temps_par_division)/points:10*temps_par_division;
    ymult=query(dso,'WFMOutpre:YMUlt?');



for iii=1:M
    disp(iii/M)
      filename=sprintf('SR360times25_500MHZ_LogP_%d.mat',iii);
    
    %%%%%%%%%%%%%MVMT BRasseur
    %clc;clear all;
serie = serial('COM1');
set(serie,'BaudRate',9600,'DataBits',8,'StopBits',1,'Parity', 'none','InputBufferSize',512,'OutputBufferSize',512,'Terminator','ETX');
set(serie,'Timeout',0.1);
serie.ReadAsyncMode = 'continuous';
STX=char(02);

%--Ouverture du port série
fopen(serie);

%--Paramétres

%angle_rotation=360/M;                              %angle de rotation en degrés.
fbrasseur=929;                                  %frequence de marche de 50 à 2000Hz.

%--------------------
% display(['--Brasseur en cours--']);
% display([' ']);
% display(['Angle de rotation (°): ' num2str(angle_rotation)]);
% display(['Fréquence de marche (Hz): ' num2str(fbrasseur)]);

%pas=angle_rotation*2275;                        %1 degré=2275 pas
pas=91; %2275/25
%--Reset

chaine=[STX '0CI:XX'];fprintf(serie,chaine);
while (serie.BytesAvailable) == 0; end ;fscanf(serie);

%--Fréquence de marche du brasseur--

chaine=[STX '0PF' num2str(fbrasseur) ':XX'];fprintf(serie,chaine); 
while (serie.BytesAvailable) == 0;end ;fscanf(serie);

%--Lecture de la position actuelle du brasseur

% chaine=[STX '0PC?:XX'];fprintf(serie,chaine); 
% while (serie.BytesAvailable) == 0; end ;out = fscanf(serie);
% e=find(out==':');
% pos_debut=str2num(out(e(1)+1:e(2)-1));
% pos_fin=pas+pos_debut;

% display(['Position actuelle du brasseur: ' num2str(pos_debut) ' pas']); 
% display(['Avance relative de: ' num2str(pas) ' pas']);

%--Avance relative de x pas

chaine=[STX '0GR' num2str(pas) ':XX'];fprintf(serie,chaine);
while (serie.BytesAvailable) == 0; end;fscanf(serie);

%--Test de fin de brassage--

%pos_cours=pos_debut;
% while (pos_cours~=pos_fin)
%     chaine=[STX '0PC?:XX'];fprintf(serie,chaine); 
%     while (serie.BytesAvailable) == 0; end ;out = fscanf(serie);
%     e=find(out==':');
%     pos_cours=str2num(out(e(1)+1:e(2)-1));
% %     display(['Position en cours du brasseur: ' num2str(pos_cours) ' pas']);
% end
% display(['Position fin du brasseur: ' num2str(pos_cours) ' pas']);
%--Fin de brassage--

%r=0.31;
%delay=0.0;%2+(r*angle_rotation);
%display(['Attente de fin de brassage et stabilisation de 15 secondes']);
%pause(delay);

%--Fermeture du port série
fclose(serie);delete(serie);clear serie;
    
    
    
    %%%%%%%%%%%%%%AQUISITION
  
    
    
    %--Enregistrement du signal CH2

    fprintf(dso,'DATa:SOURce CH2');
    fprintf(dso,'DATa:ENCDG ASCIi');
    fprintf(dso,'WFMOutpre:BIT_Nr 16');
    fprintf(dso,'DATa:STARt 1');
    stop_points=['DATa:STOP ' num2str(points) ];
    fprintf(dso,stop_points)
    fprintf(dso,'ACQuire:STATE ON');
    
   % pause
    curve = query(dso,'CURVe?');
    pause(.01);
    wave_form1 = str2num(curve)*str2num(ymult);
    wave_form1=wave_form1';
    Amplitude_crete(iii)=max(wave_form1-mean(wave_form1));
    %C(iii)=curve;
    %save(filename,'Amplitude_crete')%,'ymult')
    clear curve
    %clear wave_form1
    % figure;
    % plot(t,wave_form1);
    % title('Amplitude du signal CH1 Oscilloscope');
    % xlabel('Temps (s)');
    % ylabel('Amplitude (V)');
    % display(['N= ' num2str(points) ' points']);
    % display(['Durée= ' num2str(duree) ' secondes']);
    % display(['Amplitude maximun= ' num2str(Amplitude_crete) ' Volts']);
end
save(filename,'Amplitude_crete')
    %--Fermeture connexion AWG


fclose(dso);delete(dso);clear dso    
save tquatre.mat 't'