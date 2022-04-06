%% INSTRUMENTAÇÃO E MEDIDAS - LABORATORIO AQUISIÇÃO DE SINAIS
% Grupo 1 L32 Daniel Dinis no. 99906, João Gonçalves no. 99995, Jorge Contente no. 102143

% Dados iniciais
A = 2; % amplitude do sinal (dada em aula)
f_sinal = 200; % frequencia do sinal (dada em aula)

% Dados iniciais para a placa
Fs = 40000; % frequencia de amostragem (dada em aula)
N_amostras = 800; % no. de amostras (dada em aula)

% Resolução temporal
Ts = 1/Fs;

% Resolução espectral
F0 = Fs/N_amostras;
T0=1/F0;

% Variável no tempo e frequência
t=(0:Ts:T0-Ts)'; 
f=(0:F0:F0*(ceil(N_amostras/2)-1));

%Informacao da placa de aquisicao APAGAR
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Nbits=12;
Amax=2;
Delta=2*Amax/(2^Nbits);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Aquisição do Sinal
%d = daq("ni");
%addinput(d,"Dev2",0:1,"Voltage");
%d.Rate = Fs;
%sn_vec =read(d,N_amostras,"OutputFormat","Matrix");
%data_t1=sn_vec(:,1); %sinal na impedância desconhecida Z 
%data_t2=sn_vec(:,2); %sinal na resistência R
%Fs=sinal.Properties.SampleRate;

% Sinais de teste 
xt1 = cos(2*pi*f_sinal*t + pi/2); %sinal da impedância
xt2 = A * cos(2*pi*f_sinal*t); %sinal da resistência
data_t1=floor(xt1/Delta)*Delta+Delta/2; %sinla na impedância desconhecida Z
data_t2=floor(xt2/Delta)*Delta+Delta/2; %sinal na resistência R
%data_t1=xt1;
%data_t2=xt2;

%% Estimação da frequência
%---> estimar frequencia do sinal 1
dataf1 = abs(fft(data_t1))/N_amostras;
[M1,Posf1]=max(dataf1(1:floor(N_amostras/2),1));
media = 0;
norm = 0;
% Caso de seja possível realizar a média ponderada:
if (Posf1>3)
    for m=Posf1-3:Posf1+3
        norm = norm + dataf1(m);
        media = media +(m-1)*dataf1(m)*F0; % (m-1)*F0 é a frequencia da harmonica e dataf(m) é a sua respectiva amplitude.
                                        %(m-1) representa o nr de "subdivisoes" até à harmonica
    end
        f_estimada1 = media/norm;
else % No caso de estar proximo da origem e não dar para fazer média ponderada:
    f_estimada1 = (Posf1-1)* F0;
end
% Período estimado
T1=1/f_estimada1; 

%---> estimar frequencia do sinal 2
dataf2 = abs(fft(data_t2))/N_amostras; 
[M2,Posf2]=max(dataf2(1:floor(N_amostras/2),1));
media = 0;
norm = 0;
% Caso seja possível efetuar a média ponderada:
if (Posf2>3)
    for m=Posf2-3:Posf2+3
        norm = norm + dataf2(m);
        media = media +(m-1)*dataf2(m)*F0; % (m-1)*F0 é a frequencia da harmonica e dataf(m) é a sua amplitude.
                                        %(m-1) representa o nr de "subdivisoes" até à harmonica
    end
        f_estimada2 = media/norm;
else % No caso de estar proximo da origem e não dar para fazer média ponderada:
    f_estimada2 = (Posf2-1)* F0;
end
% Período estimado
T2=1/f_estimada2; 

%% Calculo do Navg - media do numero de amostras util quando o numero de amostras e reduzido
%---> sinal 1
nppp = Fs/f_estimada1;   	% num de pontos por periodo			    
nperiodos=floor(N_amostras/nppp);		 % num de periodos
Navg1=nperiodos*nppp;
%---> sinal 2
nppp = Fs/f_estimada2;   	% num de pontos por periodo			    
nperiodos=floor(N_amostras/nppp);		 % num de periodos
Navg2=nperiodos*nppp;


%% Calculo do Valor eficaz do:
%---> sinal 1
data_tpower=power(data_t1,2); % Vef=sqrt(mean(abs(data_t).^2))
sum_all2=sum(data_tpower);
VrmsZ=sqrt(sum_all2/Navg1);
%---> sinal 2
data_tpower=power(data_t2,2); % Vef=sqrt(mean(abs(data_t).^2))
sum_all2=sum(data_tpower);
VrmsR=sqrt(sum_all2/Navg2);

%% Diferença de fase (dif_fase)
dataf1=fft(data_t1);
dataf2=fft(data_t2);
dif_fase = angle(dataf1(Posf1)) - angle(dataf2(Posf2));
dif_fase=dif_fase*180/pi; %para ficar em radianos

%% Cálculo da impedância
% |Z| = Vz eficaz / I, I = Vr eficaz / R
% logo |Z| = Vz eficaz * R / Vr eficaz 
R=100;
abs_Z = (VrmsZ/VrmsR)*abs(R);
arg_Z = dif_fase + angle(R); %o angulo de R neste caso vai ser nulo obviamente 

%% Criar gráfico para visualização
Amplitude=max(data_t1); %amplitude real
if (Amplitude<max(data_t2))
    Amplitude=max(data_t2);
end

subplot(111);
plot(t, data_t1,'r', t, data_t2, 'b'); 
str=sprintf('Frequência do sinal de R: %g, Frequência do sinal de Z: %g, Valor eficaz de R: %g, Valor eficaz de Z: %g, Diferença de fase dos dois sinais: %g \n Número de Amostras: %g, Frequência de amostragem: %g, Alcance: [-%g, %g] V',f_estimada2, f_estimada1, VrmsR, VrmsZ, dif_fase, N_amostras, Fs, Amax, Amax);
title(str);
xlabel('t [s]')
xl = get(gca,'xlabel');
set(xl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
ylabel('Tensão [V]')
yl = get(gca,'ylabel');
set(yl,'FontName','Arial','FontSize',9,'FontWeight','bold');
legend('Sinal 1: Tensão na impedância desconhecida', 'Sinal 2: Tensão na resistência conhecida')
axis([0 T2*4 -1.1*Amplitude 1.1*Amplitude]) %[xmin xmax  ymin ymax]
