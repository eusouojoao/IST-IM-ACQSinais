%% INSTRUMENTAÇÃO E MEDIDAS - LABORATORIO AQUISIÇÃO DE SINAIS
%Grupo 1 L32 Daniel Dinis no. 99906, João Gonçalves no. 99995, Jorge Contente no. 102143

% Dados iniciais
A = 3;          %amplitude do sinal (dada em aula)
f_sinal = 2000; %frequência do sinal (dada em aula)

% Dados iniciais para a placa
N_amostras = 21000; % no. de amostras (dada em aula)
Fs = 40001;         % frequência de amostragem (dada em aula)

% Resolução temporal
Ts = 1/Fs;

% Resolução espectral
F0 = Fs/N_amostras;
T0=1/F0;

% Variável no tempo e frequência
t=(0:Ts:T0-Ts)'; 
f=(0:F0:F0*(ceil(N_amostras/2)-1));

% Informação da placa de aquisicao APAGAR
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Nbits=14;
Amax=10;
Delta=2*Amax/(2^Nbits);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% Preambulo -- Calculos dos valores pedidos
N_ciclos=1000;
f_tot=zeros(N_ciclos,1);

for l=1:N_ciclos  
    
    %% Aquisição do Sinal
    %d = daq('ni');
    %addinput(d,"Dev2","ai0","Voltage"); %daqlist
    %atribuição das especificações da sessão 
    %d.Rate = Fs;
    %inicio da aquisição
    %data_t = read(d,N_amostras,"OutputFormat","Matrix"); % Adquirir "N_amostras"  e colocar os dados na variável data_t.
    %Fs=sinal.Properties.SampleRate;

    %Funções de teste 
    xt=A*cos(2*pi*f_sinal*t+2*pi*rand(1));
    data_t=floor(xt/Delta)*Delta+Delta/2;
    %data_t=xt;
    
    % Estimação da frequência
    dataf = abs(fft(data_t))/N_amostras; 

    [M,Posf]=max(dataf(1:floor(N_amostras/2)));

    media = 0;
    norm = 0;

    %Na eventualidade da possibilidade de calcular a média ponderada:
    if (Posf>3)
        for m=Posf-3:Posf+3
            norm = norm + dataf(m);
            media = media +(m-1)*dataf(m)*F0;   
            % (m-1)*Delta_f é a frequencia da harmonica e dataf(m) é a sua amplitude.
            % (m-1) representa o nr de "subdivisoes" até à harmonica
        end
            f_tot(l) = media/norm;
    else 
        % No caso de estar proximo da origem e não dar para fazer média ponderada:
        f_tot(l) = (Posf-1)* F0;
    end

    %Calcular espetro unilateral
    dataf=fft(data_t);
    dataf_b=(abs(dataf)).^2/(N_amostras*N_amostras);    %fazer o espectro de potencia bilateral
    dataf_uni=2*dataf_b(1:floor((N_amostras+1)/2));     %fazer o espetro de potencia unilateral
    dataf_uni(1)= dataf_uni(1)/2;
    if (rem(N_amostras,2)==0)
        dataf_uni(N_amostras/2)= dataf_uni(N_amostras/2)/2;
    end
    if l == 1
        dataf_avg = dataf_uni;
    else
        dataf_avg = dataf_avg + dataf_uni;
    end
end

f_estimada=mean(f_tot);
dataf_avg = dataf_avg / N_ciclos;

% Período estimado
T=1/f_estimada; 

% Cálculo do Navg - média do numero de amostras; util quando o numero de amostras e reduzido
nppp = Fs/f_estimada;   	                % no. de pontos por periodo			   
nperiodos=floor(N_amostras/nppp);	        % num de periodos
Navg=nperiodos*nppp;

% Cálculo da média 
sum_all=sum(data_t);
media=sum_all/Navg;

% Cálculo do Valor eficaz  
data_tpower=power(data_t,2);    % Vef=sqrt(mean(abs(data_t).^2))
sum_all2=sum(data_tpower);
Vrms=sqrt(sum_all2/Navg);


%% Cálculo do ruído médio eficaz
% Naturalmente é necessário retirar as harmónicas do espectro médio de potência
% no nosso caso temos uma sinusoide, bastava retirar a fundamental, 
% mas o raciocinio está generalizado (no for loop anterior removemos todas as harmónicas do Gn)

Gn = dataf_avg;                             %salvaguardar GMed
HN = floor(N_amostras/2 - 1)*F0/f_estimada; % no. de harmónicas
for i=1:HN
    Gn(1 + round(i * f_estimada/F0)) = 0; 
end
Gm = power(10,10*log10(Gn)/10);

nQ = sum(Gm);        % ruído de quantização (ruido medio)
%nQ=nQ-Gm(Posf);     % quando temos apenas uma sinusoide (apenas um dirac no espetro de potencia unilateral)
nQ = sqrt(nQ); 


%% Cálculo de SINAD e no. bits estimado da placa de aquisição (ENOB)
SINAD=dataf_avg(Posf)/(sum(dataf_avg)-dataf_avg(Posf));
SINAD=10*log10(SINAD);
% Calculo de Numero de bits estimado da placa, 
% chamamos à variavel Nbits_Est, mas é o ENOB tecnicamente
Nbits_Est=(SINAD-1.76-10*log10(A^2/Amax^2))/6.02;

%% Criar gráficos para visualização

Amplitude=max(data_t);      % amplitude real
subplot(211);
plot(t, data_t, 'color', [0 0.5 1]); 
str=sprintf('Estimativa do ruído eficaz equivalente: %g, \n Valor da SINAD: %gdB, Frequência: %g, \n Estimativa do número efetivo de bits (ENOB): %g, \n Valor eficaz: %g, Número de Amostras: %g,\n Frequência de amostragem: %g, Alcance: [-%g, %g] V', nQ, SINAD, Nbits_Est ,f_estimada, Vrms, N_amostras, Fs, Amax, Amax);
title(str);
xlabel('t [s]')
xl = get(gca,'xlabel');
set(xl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
ylabel('Tensão [V]')
yl = get(gca,'ylabel');
set(yl,'FontName','Arial','FontSize',9,'FontWeight','bold');
legend('Gráfico temporal')
axis([0 T*5 -1.1*Amplitude 1.1*Amplitude]) %[xmin xmax  ymin ymax]

subplot(212);
max2=max(10*log10(abs(dataf_avg)))+20;          %valor max no plot
min2=1.1*min(10*log10(abs(dataf_avg)));         %valor min no plot
plot(f, 10*log10(dataf_avg), 'color',[1 0 0]);  %plot em db
%plot(f, dataf_uni, 'color',[1 0 0]);           %plot em V^2
xl = get(gca,'xlabel');
set(xl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
ylabel('Tensão [dB V]')
yl = get(gca,'ylabel');
set(yl,'FontName','Arial','FontSize',9,'FontWeight','bold');   
legend('Espectro')
axis([0 Fs/2  min2 max2])
