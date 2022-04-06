%% INSTRUMENTAÇÃO E MEDIDAS - LABORATORIO AQUISIÇÃO DE SINAIS
% Grupo 1 L32 Daniel Dinis no. 99906, João Gonçalves no. 99995, Jorge Contente no. 102143

% Dados iniciais
A = 3.5; %amplitude do sinal (dada em aula)
f_sinal = 900; %frequência do sinal (dada em aula)
N_amostras = 5000; %no. de amostras (dada em aula)
Fs = 36000; %frequência de amostragem (dada em aula)

% Informacao da placa de aquisicao (APAGAR)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Nbits=12;
Amax=10;
Delta=2*Amax/(2^Nbits); %dado para testes
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Aquisição do Sinal
%d = daq('ni');
%addinput(d,"Dev2","ai0","Voltage");
%atribuição das especificações da sessão
%d.Rate = Fs;
%inicio da aquisição
%data_t = read(d,N_amostras,"OutputFormat","Matrix"); % Adquirir "N_amostras"  e colocar os dados na variável data_t.
%Fs=sinal.Properties.SampleRate;


%% Estimar a frequência do sinal
% Resolucao temporal
Ts = 1/Fs;

% Resolucao espectral
F0 = Fs/N_amostras;
T0=1/F0;

% Variável no tempo e frequência
t=(0:Ts:T0-Ts)'; 
f=(0:F0:F0*(ceil(N_amostras/2)-1));


% Funções de teste
%xt=A*cos(2*pi*f_sinal*t);
xt=A*sawtooth(2*pi*f_sinal*t,0.5);
%xt=A*square(2*pi*f_sinal*t);

% Simular placa de aquisição
data_t=floor(xt/Delta)*Delta+Delta/2;
%data_t=xt;


% Obter máximo e posição da fft
dataf = abs(fft(data_t))/N_amostras; 

[M,Posf]=max(dataf(1:floor(N_amostras/2),1));

media = 0;
norm = 0;

% Fazer média ponderada:
if (Posf>3)
    for m=Posf-3:Posf+3
        norm = norm + dataf(m);
        media = media +(m-1)*dataf(m)*F0; % (m-1)*F0 é a frequencia da harmonica e dataf(m) é a sua respectiva amplitude.
                                        %(m-1) representa o no. de "subdivisoes" até à harmónica
    end
        f_estimada = media/norm;
else % No caso de estar proximo da origem e não dar para fazer média ponderada:
    f_estimada = (Posf-1)* F0;
end

% Periodo estimado
T=1/f_estimada; 

%% Calculo da media
% Calculo do Navg - media do numero de amostras (reduz o espalhamento)

nppp = Fs/f_estimada;   	% num de pontos por periodo			    
nperiodos=floor(N_amostras/nppp);		 % num de periodos
Navg=nperiodos*nppp;

% Media
sum_all=sum(data_t);
media=sum_all/Navg;

%% Calculo do Valor eficaz  

data_tpower=power(data_t,2); % Vrms=sqrt(mean(abs(data_t).^2))
sum_all2=sum(data_tpower);
Vrms=sqrt(sum_all2/Navg);

%% Calculo das harmónicas
% Descobrir o número total de harmónicas
PosfH=f_estimada/F0; %Posf da fundamental +1
n_harmonicas=floor((N_amostras/2-2)/PosfH); %número de harmónicas
% Vetores para armanezar o resultado final
amplitude_harmonica=zeros(1, n_harmonicas+1);
amplitude_harmonica_dB=zeros(1, n_harmonicas+1);

% Espetro do sinal, normalizado e unilaterlizado
dataf = abs(fft(data_t))/N_amostras;
dataf=2*dataf(1:floor(N_amostras/2)); 
dataf(1)=dataf(1)/2;

% Colocar harmónicas nos vetores
for l=1:n_harmonicas+1
    amplitude_harmonica(l)=dataf(round((l-1)*PosfH)+1)/sqrt(2); %divide-se por sqrt(2) para ficar valor eficaz
    amplitude_harmonica_dB(l)=20*log10(dataf(round((l-1)*PosfH)+1));
end

% Nota: 1a. componente do vetor contém a compontente DC, 
% a 1a. harmónica está na posição 2 do vetor
for n=2:n_harmonicas+1 %print amplitude harmonica
    %fprintf(1,'Amplitude da harmonica de tensao nr %d = %.4f \n',n-1,amplitude_harmonica(n));
end
%fprintf('\n');

for n=2:n_harmonicas+1 %print amplitude harmonica dB
    %fprintf(1,'Amplitude em dB da harmonica de tensao nr %d = %.4f \n',n-1,amplitude_harmonica_dB(n));
end

%% Calculo da distorção harmónica total (THD)
sum_minor=sum(amplitude_harmonica(3:n_harmonicas+1).^2);
THD = 20*log10(sqrt(sum_minor/amplitude_harmonica(2)^2)); %pois a amplitude_harmonica(1) tem a componente DC 

%% Criar gráfico temporal
Amplitude=max(data_t); %amplitude real
subplot(211);
plot(t, data_t, 'color', [0 0.5 1]); 

str=sprintf('Frequência: %g, Valor médio: %g, Valor eficaz: %g,\n Número de Amostras: %g, Frequência de amostragem: %g, Alcance: [%g, %g] V', f_estimada, media, Vrms, N_amostras, Fs, -Amax, Amax);
title(str);
xlabel('t [s]')
xl = get(gca,'xlabel');
set(xl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
ylabel('Tensão [V]')
yl = get(gca,'ylabel');
set(yl,'FontName','Arial','FontSize',9,'FontWeight','bold');
legend('Gráfico temporal')
axis([0 T*5 -1.1*Amplitude 1.1*Amplitude]) %[xmin xmax  ymin ymax]
                                           %T*5 para serem 5 periodos

%% Gráfico de potência unilateral para visualização
subplot(212);
% Calcular espetro de potência unilateral
dataf=fft(data_t);
dataf_b=(abs(dataf)).^2/(N_amostras*N_amostras); %fazer o espectro de potência bilateral
dataf_uni=2*dataf_b(1:floor((N_amostras+1)/2)); %fazer o espetro de potência unilateral
dataf_uni(1)= dataf_uni(1)/2; %dataf_uni = dataf_b, se m == 0
if (rem(N_amostras,2)==0) %dataf_uni = dataf_b, se m == N_amostras for par
    dataf_uni(N_amostras/2)= dataf_uni(N_amostras/2)/2;
end;

max2=max(10*log10(abs(dataf_uni)))+10; %valor max no plot
min2=1.1*min(10*log10(abs(dataf_uni))); %valor min no plot

plot(f, 10*log10(dataf_uni),'color',[0 0.5 0.2]); %plot em db
%plot(f, dataf_uni, 'color',[0 0.5 0.2]); %plot em V^2

xlabel('f [Hz]')
xl = get(gca,'xlabel');
set(xl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
ylabel('Tensão [dB V]')
yl = get(gca,'ylabel');
set(yl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
legend('Espectro')
axis([0 Fs/2  min2 max2])
