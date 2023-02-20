close all; clear; clc;
%Definisanje osnovnih postavki simulacije
Fs=384; %Frekvencija odabiranja kojom se generisu signali
Fg=6; %Frekvencija signala
Tsample=1/Fs; %Perioda odabiranja
Nsegm=16; #Ukupan broj segmenata
SizeSegm=4096; %Duzina segmenata
TrajanjeSegmenta=SizeSegm*Tsample; %TrajanjeSegmenta u sekundama
Nsamples=Nsegm*SizeSegm; %Ukupan broj odbiraka za svaki od generisanih signala

%Vremenska osa
tosa=0:Tsample:(Nsamples-1)*Tsample;

%Generisanje periodicnih signala
Nsig=4; %Broj komponenti slozenoperiodicnog signala
A=[2 1 0.5 0.25]; #Amplitude komponenti slozenoperiodicnog signala
fpps=[1.5 3 4.5 6];
Ssps(1:Nsamples,1) = zeros(Nsamples,1);
    Ssps(1:Nsamples,1) = Ssps(1:Nsamples,1) + 2;    % Dodavanje DC komponente
    for i = 1 : Nsig                              % Dodavanje prostoperiodicnih komponenti
        Ssps(1:Nsamples,1) = Ssps(1:Nsamples,1) + A(i)*cos(2*pi*fpps(i)*tosa)';
    end;


 % Crtanje vremeskog oblika svih generisanih signala
    figure(1) % Signal 1
    hold on; grid on; box on;
    subplot(2,1,1)
    prikaz = 0;
    prikaz = Ssps(:,1);
    plot(tosa,prikaz,'b')
    title('Slozenoperiodican signal')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 5 min(prikaz)-0.1 max(prikaz)+0.5])




    % Definisanje modulišuceg signala za AM modulaciju

    um(1:Nsamples,1) = Ssps(1:Nsamples,1);


    % Definisanje ucestanosti mosioca i signala nosioca na predaji
    fcarrier = Fs/4;        % Ucestanost nosioca 0.25 x frekvencija odabiranja
    Icarrier(1:Nsamples,1) = cos(2*pi*fcarrier*tosa);   % Nosioca u fazi
    Qcarrier(1:Nsamples,1) = sin(2*pi*fcarrier*tosa);   % Nosioc u kvadraturi

    % Definisanje ucestanosti nosioca na prijemu (za koherentnu demodulaciju) - Idealna sinhronizacija faze
    LOIcarrier(1:Nsamples,1) = 2*cos(2*pi*fcarrier*tosa);   % Lokalno generisani nosioc u fazi na prijemu - Idealna sinhronizacija faze
    LOQcarrier(1:Nsamples,1) = 2*sin(2*pi*fcarrier*tosa);   % Lokalno generisani nosioc u kvadraturi na prijemu - Idealna sinhronizacija faze

    % Definisanje ucestanosti mosioca na prijemu (za koherentnu demodulaciju) - Neidealna sinhronizacija faze
    teta1 = pi/36; teta2 = pi/18; teta3 = pi/2;
    LOIcarrierPO1(1:Nsamples,1) = 2*cos(2*pi*fcarrier*tosa + teta1);   % Lokalno generisani nosioc u fazi na prijemu - greska sinhron. faze teta1
    LOQcarrierPO1(1:Nsamples,1) = 2*sin(2*pi*fcarrier*tosa + teta1);   % Lokalno generisani nosioc u kvadraturi na prijemu - greska sinhron. faze teta1

    LOIcarrierPO2(1:Nsamples,1) = 2*cos(2*pi*fcarrier*tosa + teta2);   % Lokalno generisani nosioc u fazi na prijemu - greska sinhron. faze teta2
    LOQcarrierPO2(1:Nsamples,1) = 2*sin(2*pi*fcarrier*tosa + teta2);   % Lokalno generisani nosioc u kvadraturi na prijemu - greska sinhron. faze teta2

    LOIcarrierPO3(1:Nsamples,1) = 2*cos(2*pi*fcarrier*tosa + teta3);   % Lokalno generisani nosioc u fazi na prijemu - greska sinhron. faze teta3
    LOQcarrierPO3(1:Nsamples,1) = 2*sin(2*pi*fcarrier*tosa + teta3);   % Lokalno generisani nosioc u kvadraturi na prijemu - greska sinhron. faze teta3

    % Definisanje ucestanosti mosioca na prijemu (za koherentnu demodulaciju) - greška sinhronizacije ucestanosti (frekvencijski offset)
    freqoffset = 5;         % Greška ucestanosti 5Hz
    LOIcarrierFO(1:Nsamples,1) = 2*cos(2*pi*(fcarrier + freqoffset)*tosa);  % Lokalno generisani nosioc u fazi na prijemu - greska sinhron. ucestanosti
    LOQcarrierFO(1:Nsamples,1) = 2*sin(2*pi*(fcarrier + freqoffset)*tosa);  % Lokalno generisani nosioc u kvadraturi na prijemu - greska sinhron. ucestanosti
    % Generisanje LPF granicne ucestanosti Fs reda 256
    LPForder1=256;
    BPForder=256;
    LPF = fir1(LPForder1,Fg/(Fs/2));
    % Generisanje BPF filtara reda 256
    BPF1 = fir1(BPForder,[Fg/(Fs/2) 2*Fg/(Fs/2)]);
    BPF2 = fir1(BPForder,[Fg/Fs 2*Fg/Fs]);


    % CAM signal - koherentna i kvadraturna demodulacija
    Apps=5.75; %Maksimalna vrednost-Kada su kosinusi 1
    m0 =[1/2 1/4 5/4];    % Stepen modulacije
    SNR=20;  %Zadati SNR za AWGN funkciju
    CAM_Rx=zeros(1:Nsamples,6);
    for i=1:3
        CAM_Rx(1:Nsamples,i) = (1+m0(i)/Apps*um(1:Nsamples,1)).*Icarrier(1:Nsamples,1);   % CAM modulisani signala na prijemu
        CAM_Rx(1:Nsamples,i)=awgn(CAM_Rx(1:Nsamples,i),SNR,'measured');    % Dodavanje ABGŠ sa zadatim odnosom SNR
        % Koherentna demodulacija CAM - sa idealnom sinhrnonizacijom i neidealnom sinhronizacijom - signal nakon množenja sa lokalno generisanim nosiocem
        CAM_DemI(1:Nsamples,i) = CAM_Rx(1:Nsamples,i).*LOIcarrier(1:Nsamples,1);        % Idealna sinhronizacija
        CAM_DemIPO1(1:Nsamples,i) = CAM_Rx(1:Nsamples,i).*LOIcarrierPO1(1:Nsamples,1);  % Neidealna sinhronizacija - teta1
        CAM_DemIPO2(1:Nsamples,i) = CAM_Rx(1:Nsamples,i).*LOIcarrierPO2(1:Nsamples,1);  % Neidealna sinhronizacija - teta2
        CAM_DemIPO3(1:Nsamples,i) = CAM_Rx(1:Nsamples,i).*LOIcarrierPO3(1:Nsamples,1);  % Neidealna sinhronizacija - teta3

        % Koherentna demodulacija CAM - signal na izlazu prijemnika (filtriran LPF filtrom fg)
        CAM_izlazI(:,i) = filter(LPF,1,CAM_DemI(1:Nsamples,i));          % Idealna sinhronizacija
        CAM_izlazIPO1(:,i) = filter(LPF,1,CAM_DemIPO1(1:Nsamples,i));    % Neidealna sinhronizacija - teta1
        CAM_izlazIPO2(:,i) = filter(LPF,1,CAM_DemIPO2(1:Nsamples,i));    % Neidealna sinhronizacija - teta2
        CAM_izlazIPO3(:,i) = filter(LPF,1,CAM_DemIPO3(1:Nsamples,i));    % Neidealna sinhronizacija - teta3

        % Odbacivanje jednosmerne komponente za koherenetan prijemnik
        CAM_izlazIac(:,i) = filter(BPF1,1,CAM_DemI(1:Nsamples,1));          % Idealna sinhronizacija
        CAM_izlazIPO1ac(:,i) = filter(BPF1,1,CAM_DemIPO1(1:Nsamples,1));    % Neidealna sinhronizacija - teta1
        CAM_izlazIPO2ac(:,i) = filter(BPF1,1,CAM_DemIPO2(1:Nsamples,1));    % Neidealna sinhronizacija - teta2
        CAM_izlazIPO3ac(:,i) = filter(BPF1,1,CAM_DemIPO3(1:Nsamples,1));    % Neidealna sinhronizacija - teta3
        CAM_izlazIac(:,i+3) = filter(BPF2,1,CAM_DemI(1:Nsamples,1));          % Idealna sinhronizacija
        CAM_izlazIPO1ac(:,i+3) = filter(BPF2,1,CAM_DemIPO1(1:Nsamples,1));    % Neidealna sinhronizacija - teta1
        CAM_izlazIPO2ac(:,i+3) = filter(BPF2,1,CAM_DemIPO2(1:Nsamples,1));    % Neidealna sinhronizacija - teta2
        CAM_izlazIPO3ac(:,i+3) = filter(BPF2,1,CAM_DemIPO3(1:Nsamples,1));    % Neidealna sinhronizacija - teta3

    endfor
    ;
    % Koherentna demodulacija CAM - sa idealnom sinhrnonizacijom i neidealnom sinhronizacijom - signal nakon množenja sa lokalno generisanim nosiocem
    CAM_DemI(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOIcarrier(1:Nsamples,1);        % Idealna sinhronizacija
    CAM_DemIPO1(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOIcarrierPO1(1:Nsamples,1);  % Neidealna sinhronizacija - teta1
    CAM_DemIPO2(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOIcarrierPO2(1:Nsamples,1);  % Neidealna sinhronizacija - teta2
    CAM_DemIPO3(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOIcarrierPO3(1:Nsamples,1);  % Neidealna sinhronizacija - teta3

    % Koherentna demodulacija CAM - signal na izlazu prijemnika (filtriran LPF filtrom fg = 12Hz)
    CAM_izlazI(:,1) = filter(LPF,1,CAM_DemI(1:Nsamples,1));          % Idealna sinhronizacija
    CAM_izlazIPO1(:,1) = filter(LPF,1,CAM_DemIPO1(1:Nsamples,1));    % Neidealna sinhronizacija - teta1
    CAM_izlazIPO2(:,1) = filter(LPF,1,CAM_DemIPO2(1:Nsamples,1));    % Neidealna sinhronizacija - teta2
    CAM_izlazIPO3(:,1) = filter(LPF,1,CAM_DemIPO3(1:Nsamples,1));    % Neidealna sinhronizacija - teta3

    % Odbacivanje jednosmerne komponente za koherenetan prijemnik
    CAM_izlazIac(:,1) = filter(BPF1,1,CAM_DemI(1:Nsamples,1));          % Idealna sinhronizacija
    CAM_izlazIPO1ac(:,1) = filter(BPF1,1,CAM_DemIPO1(1:Nsamples,1));    % Neidealna sinhronizacija - teta1
    CAM_izlazIPO2ac(:,1) = filter(BPF1,1,CAM_DemIPO2(1:Nsamples,1));    % Neidealna sinhronizacija - teta2
    CAM_izlazIPO3ac(:,1) = filter(BPF1,1,CAM_DemIPO3(1:Nsamples,1));    % Neidealna sinhronizacija - teta3


    % Kvadraturna (nekoherentna) demodulacija CAM - grana u kvadraturi
    CAM_DemQ(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOQcarrier(1:Nsamples,1);        % Idealna sinhronizacija
    CAM_DemQPO1(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOQcarrierPO1(1:Nsamples,1);  % Neidealna sinhronizacija - teta1
    CAM_DemQPO2(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOQcarrierPO2(1:Nsamples,1);  % Neidealna sinhronizacija - teta2
    CAM_DemQPO3(1:Nsamples,1) = CAM_Rx(1:Nsamples,1).*LOQcarrierPO3(1:Nsamples,1);  % Neidealna sinhronizacija - teta3

    CAM_izlazQ(:,1) = filter(LPF,1,CAM_DemQ(1:Nsamples,1));          % Idealna sinhronizacija
    CAM_izlazQPO1(:,1) = filter(LPF,1,CAM_DemQPO1(1:Nsamples,1));    % Neidealna sinhronizacija - teta1
    CAM_izlazQPO2(:,1) = filter(LPF,1,CAM_DemQPO2(1:Nsamples,1));    % Neidealna sinhronizacija - teta2
    CAM_izlazQPO3(:,1) = filter(LPF,1,CAM_DemQPO3(1:Nsamples,1));    % Neidealna sinhronizacija - teta3

    % Izlaz kvadraturnog demodulatora
    CAM_izlazQAM = sqrt(CAM_izlazI.^2 + CAM_izlazQ.^2);              % Idealna sinhronizacija
    CAM_izlazQAMPO1 = sqrt(CAM_izlazIPO1.^2 + CAM_izlazQPO1.^2);     % Neidealna sinhronizacija - teta1
    CAM_izlazQAMPO2 = sqrt(CAM_izlazIPO2.^2 + CAM_izlazQPO2.^2);     % Neidealna sinhronizacija - teta2
    CAM_izlazQAMPO3 = sqrt(CAM_izlazIPO3.^2 + CAM_izlazQPO3.^2);     % Neidealna sinhronizacija - teta3

    % Odbacivanje jednosmerne komponente za kvadraturni prijemnik
    for i=1:3
      CAM_izlazQAMac(:,i) = filter(BPF1,1,CAM_izlazQAM(1:Nsamples,i));          % Idealna sinhronizacija
      CAM_izlazQAMPO1ac(:,i) = filter(BPF1,1,CAM_izlazQAMPO1(1:Nsamples,i));    % Neidealna sinhronizacija - teta1
      CAM_izlazQAMPO2ac(:,i) = filter(BPF1,1,CAM_izlazQAMPO2(1:Nsamples,i));    % Neidealna sinhronizacija - teta2
      CAM_izlazQAMPO3ac(:,i) = filter(BPF1,1,CAM_izlazQAMPO3(1:Nsamples,i));    % Neidealna sinhronizacija - teta3
      CAM_izlazQAMac(:,i+3) = filter(BPF2,1,CAM_izlazQAM(1:Nsamples,i));          % Idealna sinhronizacija
      CAM_izlazQAMPO1ac(:,i+3) = filter(BPF2,1,CAM_izlazQAMPO1(1:Nsamples,i));    % Neidealna sinhronizacija - teta1
      CAM_izlazQAMPO2ac(:,i+3) = filter(BPF2,1,CAM_izlazQAMPO2(1:Nsamples,i));    % Neidealna sinhronizacija - teta2
      CAM_izlazQAMPO3ac(:,i+3) = filter(BPF2,1,CAM_izlazQAMPO3(1:Nsamples,i));    % Neidealna sinhronizacija - teta3
    endfor
    ;



   % Petlja za procenu spektra signala
    CAM_Tx=CAM_Rx; %Nema slabljenja na liniji
    prompt = "What is the original value? ";
    Nfft3 = input(prompt); %Uneti broj tacaka [1024 2048 4096]
    Nsegfft3 = input(prompt); %Uneti broj segmenata [64 32 16]
    fosa4 = -Fs/2:Fs./Nfft3:Fs/2 - Fs./Nfft3;

    signal = 0; signalTx = 0; signalRx = 0;
    signalDemI = 0; signalILPF = 0;
    signalDemFO = 0; signalFOLPF = 0;
    signal(1:Nsamples,1) = um(1:Nsamples,1);
    signalTx(1:Nsamples,1) = CAM_Tx(1:Nsamples,1);
    signalRx(1:Nsamples,1) = CAM_Rx(1:Nsamples,1);
    signalDemI(1:Nsamples,1) = CAM_DemI(1:Nsamples,1);
    signalILPF(1:Nsamples,1) = CAM_izlazIac(1:Nsamples,1);

    AmpSpekCAM(1:Nfft3,1:5) = zeros(Nfft3,5);
    for kkk = 1:Nsegfft3
        trensig = 0; trenspekt = 0;
        trensig1 = 0; trenspekt1 = 0;
        trensig2 = 0; trenspekt2 = 0;
        trensig3 = 0; trenspekt3 = 0;
        trensig4 = 0; trenspekt4 = 0;

        trensig(1:Nfft3,1) = signal((kkk-1)*Nfft3+1:kkk*Nfft3,1);
        trenspekt(1:Nfft3,1) = fft(trensig(1:Nfft3,1),Nfft3)/Nfft3;

        trensig1(1:Nfft3,1) = signalTx((kkk-1)*Nfft3+1:kkk*Nfft3,1);
        trenspekt1(1:Nfft3,1) = fft(trensig1(1:Nfft3,1),Nfft3)/Nfft3;

        trensig2(1:Nfft3,1) = signalRx((kkk-1)*Nfft3+1:kkk*Nfft3,1);
        trenspekt2(1:Nfft3,1) = fft(trensig2(1:Nfft3,1),Nfft3)/Nfft3;

        trensig3(1:Nfft3,1) = signalDemI((kkk-1)*Nfft3+1:kkk*Nfft3,1);
        trenspekt3(1:Nfft3,1) = fft(trensig3(1:Nfft3,1),Nfft3)/Nfft3;

        trensig4(1:Nfft3,1) = signalILPF((kkk-1)*Nfft3+1:kkk*Nfft3,1);
        trenspekt4(1:Nfft3,1) = fft(trensig4(1:Nfft3,1),Nfft3)/Nfft3;

        AmpSpekCAM(1:Nfft3,1) = AmpSpekCAM(1:Nfft3,1) + (abs(trenspekt(1:Nfft3,1)).^2);
        AmpSpekCAM(1:Nfft3,2) = AmpSpekCAM(1:Nfft3,2) + (abs(trenspekt1(1:Nfft3,1)).^2);
        AmpSpekCAM(1:Nfft3,3) = AmpSpekCAM(1:Nfft3,3) + (abs(trenspekt2(1:Nfft3,1)).^2);
        AmpSpekCAM(1:Nfft3,4) = AmpSpekCAM(1:Nfft3,4) + (abs(trenspekt3(1:Nfft3,1)).^2);
        AmpSpekCAM(1:Nfft3,5) = AmpSpekCAM(1:Nfft3,5) + (abs(trenspekt4(1:Nfft3,1)).^2);
    end;
    AmpSpekCAM(1:Nfft3,1:5) = sqrt(AmpSpekCAM(1:Nfft3,1:5)/Nsegfft3/Fs*Nfft3);

    clear signal signalTx signalRx signalDemI signalILPF signalDemFO signalFOLPF
    clear trensig trenspekt trensig1 trenspekt1 trensig2 trenspekt2
    clear trensig3 trenspekt3 trensig4 trenspekt4
    clear trensig5 trenspekt5 trensig6 trenspekt6

    % Prikaz - Vremenski oblici signala za idealnu sinhronizaciju prijemnika

    % Crtanje vremenskog oblika signala u razlicitim tackama sistema
    figure(2) % CAM signali - Idealna sinhronizacija

    subplot(3,1,1)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = um(1:Nsamples,1);
    plot(tosa,prikaz,'b')
    title('CAM modulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 3 min(prikaz)-0.1 max(prikaz)+0.5])

    prikaz = 0;
    prikaz = CAM_Tx(1:Nsamples,1);
    plot(tosa,prikaz,'m')

    prikaz = 0;
    prikaz = CAM_Rx(1:Nsamples,1);
    plot(tosa,prikaz,'r')

    legend('Modulisuci signal za m0=1/2','AM-DSB signal - Izlaz predajnika','AM-DSB signal - Ulaz prijemnika')

    subplot(3,1,2)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = um(1:Nsamples,1);
    plot(tosa,prikaz,'b')
    title('CAM modulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 3 min(prikaz)-0.1 max(prikaz)+0.5])

    prikaz = 0;
    prikaz = CAM_Tx(1:Nsamples,2);
    plot(tosa,prikaz,'m')

    prikaz = 0;
    prikaz = CAM_Rx(1:Nsamples,2);
    plot(tosa,prikaz,'r')

    legend('Modulisuci signal za m0=1/4','AM-DSB signal - Izlaz predajnika','AM-DSB signal - Ulaz prijemnika')

    subplot(3,1,3)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = um(1:Nsamples,1);
    plot(tosa,prikaz,'b')
    title('CAM modulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 3 min(prikaz)-0.1 max(prikaz)+0.5])

    prikaz = 0;
    prikaz = CAM_Tx(1:Nsamples,3);
    plot(tosa,prikaz,'m')

    prikaz = 0;
    prikaz = CAM_Rx(1:Nsamples,3);
    plot(tosa,prikaz,'r')

    legend('Modulisuci signal za m0=5/4','AM-DSB signal - Izlaz predajnika','AM-DSB signal - Ulaz prijemnika')

    clear prikaz
   %%Demodulisan
    figure(3)
    subplot(4,1,1)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = CAM_izlazIac(1:Nsamples,1);
    plot(tosa,prikaz,'r')
    title('CAM demodulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 5 min(prikaz)-0.1 max(prikaz)+0.1])



    legend('Demodulisan signal na izlazu')

    subplot(4,1,2)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = CAM_izlazIac(1:Nsamples,2);
    plot(tosa,prikaz,'r')
    title('CAM demodulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 5 min(prikaz)-0.1 max(prikaz)+0.1])

    legend('Demodulisan signal na izlazu')

    subplot(4,1,3)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = CAM_izlazIac(1:Nsamples,3);
    plot(tosa,prikaz,'r')
    legend('Demodulisan signal na izlazu')

    title('CAM demodulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 5 min(prikaz)-0.1 max(prikaz)+0.1])


    subplot(4,1,4)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = CAM_izlazIac(1:Nsamples,4);
    plot(tosa,prikaz,'r')
    title('CAM demodulacija')
    xlabel('Vreme [s]')
    ylabel('Vremenski oblik signala')
    axis([0 5 min(prikaz)-0.1 max(prikaz)+0.1])

    legend('Demodulisan signal na izlazu')

    clear prikaz

    % Crtanje amp. spektara signala u razlicitim tackama sistema
    figure(4)
    hold on

    subplot(5,1,1)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = 20*log10(AmpSpekCAM(1:Nfft3,1));
    plot(fosa4,fftshift(prikaz),'r')
    title('Modulisuci signal')
    xlabel('Ucestanost [Hz]')
    ylabel('Amp. spektar signala')
    axis([-600 600 min(prikaz)-0.1 max(prikaz)+0.5])

    subplot(5,1,2)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = 20*log10(AmpSpekCAM(1:Nfft3,2));
    plot(fosa4,fftshift(prikaz),'r')
    title('CAM signal - Izlaz predajnika')
    xlabel('Ucestanost [Hz]')
    ylabel('Amp. spektar signala')
    axis([-600 600 min(prikaz)-0.1 max(prikaz)+0.5])

    subplot(5,1,3)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = 20*log10(AmpSpekCAM(1:Nfft3,3));
    plot(fosa4,fftshift(prikaz),'r')
    title('CAM signal - Ulaz prijemnika')
    xlabel('Ucestanost [Hz]')
    ylabel('Amp. spektar signala')
    axis([-600 600 min(prikaz)-0.1 max(prikaz)+0.5])

    subplot(5,1,4)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = 20*log10(AmpSpekCAM(1:Nfft3,4));
    plot(fosa4,fftshift(prikaz),'r')
    title('CAM signal - Nakon mnozenja lokalnim nosiocem')
    xlabel('Ucestanost [Hz]')
    ylabel('Amp. spektar signala')
    axis([-600 600 min(prikaz)-0.1 max(prikaz)+0.5])

    subplot(5,1,5)
    hold on; grid on; box on;
    prikaz = 0;
    prikaz = 20*log10(AmpSpekCAM(1:Nfft3,5));
    plot(fosa4,fftshift(prikaz),'r')
    title('CAM signal - na izlazu prijemnika')
    xlabel('Ucestanost [Hz]')
    ylabel('Amp. spektar signala')
    axis([-600 600 min(prikaz)-0.1 max(prikaz)+0.5])

    clear prikaz

    % Prikaz - Vremenski oblici signala za neidealnu sinhronizaciju faze prijemnika
    for i=1:6
      figure(4+i) % CAM - Idealna i neidealna sinhronizacija faze
      hold on; grid on; box on;

      prikaz = 0;
      prikaz = CAM_izlazIac(1:Nsamples,i);
      plot(tosa,prikaz,'b')
      title('Demodulisan signal - Idealna i neidealna sinhronizacija faze')
      xlabel('Vreme [s]')
      ylabel('Vremenski oblik signala')
      axis([0 4 min(prikaz)-0.1 max(prikaz)+0.5])

      prikaz = 0;
      prikaz = CAM_izlazIPO1ac(1:Nsamples,i);
      plot(tosa,prikaz,'r')

      prikaz = 0;
      prikaz = CAM_izlazIPO2ac(1:Nsamples,i);
      plot(tosa,prikaz,'m')

      prikaz = 0;
      prikaz = CAM_izlazIPO3ac(1:Nsamples,i);
      plot(tosa,prikaz,'c')
    endfor



    % Vremenski oblik signala za kvadraturnu demodulaciju - Idealna i neidealna sinhronizacija faze

    % Crtanje vremenskog oblika signala u razlicitim tackama sistema
    % Crtanje vremenskog oblika za mo=[0.5 0.25 1.25] aproksimativno
for i=1:6
  figure(10+i)
      hold on; grid on; box on;
      figure(10+i)
      prikaz = 0;
      prikaz = 1.01*CAM_izlazQAMac(1:Nsamples,i);
      plot(tosa,prikaz,'b')
      title('Demodulisan signal - Idealna i neidealna sinhronizacija faze')
      xlabel('Vreme [s]')
      ylabel('Vremenski oblik signala')
      axis([0 7.5 min(prikaz)-0.1 max(prikaz)+0.5])
      figure(10+i)
      prikaz = 0;
      prikaz = CAM_izlazQAMPO1ac(1:Nsamples,i);
      plot(tosa,1.01*prikaz,'r')
      figure(10+i)
      prikaz = 0;
      prikaz = CAM_izlazQAMPO2ac(1:Nsamples,i);
      plot(tosa,0.98*prikaz,'b')
      figure(10+i)
      prikaz = 0;
      prikaz = CAM_izlazQAMPO3ac(1:Nsamples,i);
      plot(tosa,1.02*prikaz,'c')

      hold on

      legend('Demodulisan signal - Idealna sinhr. faze','Kvadraturno demodulisan signal - Greska sinhr. faze pi/36','Kvadraturno demodulisan signal - Greska sinhr. faze pi/18','Kvadraturno demodulisan signal - Greska sinhr. faze pi/2')
endfor;


