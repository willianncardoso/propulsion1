% Turbofan ideal 
% aqui corremos em pi_c (raz?o de press?o de compressor)

%zerando------------------------------------------------------------------%
clear variables
close all
clc
format long
%-------------------------------------------------------------------------%

% Dados iniciais conhecidos-----------------------------------------------%
T0   = 216.7;                   %temperatura inicial    |K|
y    = 1.4;                     %                       |J/kgK|
cp   = 1004;                    %                       |kJ/kg| 
hpr  = 42800;                              
Tt4  = 1560;                    %temperaturas no fim da camara de combustao |K|
pi_f = 2;                
M0   = 0.9;                     %M0 
pi_c = 2:2:40;                  %varredura em pi_c
alfa = [4;6;8;10;12];  %CONJUNTO FECHADO razao de bypass  
                                %(alfalinha - razao otima)

j = length(alfa);               %cada linha do grafico CONJUNTO FECHADO
n = length(pi_c);               %onde corremos
%-------------------------------------------------------------------------%

%saida--------------------------------------------------------------------%
%       empuxo = |N/(kg/s)|
%       f      = |adimensional| raz?o combust?vel/ air
%       s      = |(mg/s)/N| consumo espef?cico de combust?vel por empuxo
%       nT     = |adimensional| eficiencia termica
%       nP     = |adimensional| eficiencia propulsiva
%       n0     = |adimensional| eficiencia total
%       FR     = raz?o do empuxo espec?fico por unidade de fluxo de massa 
%                do escoamento no n?cleo com aquele do FAN
%-------------------------------------------------------------------------%

%pre-alocando matrizes com zeros para preenche-las------------------------%
tal_r   = zeros(j,n); 
tal_lambida = zeros(j,n); 
raz_1   = zeros(j,n);
raz_2   = zeros(j,n);
ftotal  = zeros(j,n);
s       = zeros(j,n);
nP      = zeros(j,n);
nT      = zeros(j,n);
n0      = zeros(j,n);
%-------------------------------------------------------------------------%



%loop de varredura em pi_c--------------------------------------------------%
for i = 1:j   
    for k = 1:n

    R           = ((y-1)/y)*cp;
    a0          = sqrt(y*R*T0);
    tal_r       = 1+((y-1)/2).*(M0.^2);
    tal_lambida = (Tt4)/(T0);
    tal_c       = (pi_c).^((y-1)/y);
    tal_f       = (pi_f).^((y-1)/y);      
    %raz_1 = V9/a0
    raz_1       = sqrt((2/(y-1)).*(tal_lambida-tal_r.*((tal_c-1)+alfa.*(tal_f-1))-((tal_lambida)./(tal_r.*tal_c))));
    %raz_3 = V19/a0
    raz_3       = sqrt((2/(y-1)).*(tal_r.*tal_f-1));
    %raz_2 = F/m0ponto EMPUXO
    raz_2       = (a0./(alfa+1)).*((raz_1-M0)+alfa.*(raz_3-M0));
    f           = ((cp*T0)/hpr).*(tal_lambida-tal_r.*tal_c)/1000;
    s           = (f./((1+alfa).*raz_2)).*1000000;
    nT          = (1-1./(tal_r.*tal_c));
    nP          = (2.*M0).*((raz_1-M0+alfa.*(raz_3-M0))./(raz_1.^2-M0.^2+alfa.*(raz_3.^2-M0.^2)));
    n0          = (nT.*nP);
    FR          = (raz_1-M0)./(raz_3-M0);
    
    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_c,raz_2)
        xlabel('pi_c')
        ylabel('empuxo')
        grid
        legend('alfa_1','alfa_2','alfa_3','alfa_4','alfa_5')
        
        subplot(2,2,2)
        plot(pi_c,f)
        xlabel('pi_c')
        ylabel('f')
        grid
  
        subplot(2,2,3)
        plot(pi_c,s)
        xlabel('pi_c')
        ylabel('s')
        grid
        
        subplot(2,2,4)
        plot(pi_c,nT)
        xlabel('pi_c')
        ylabel('nT')
        grid
        
        figure(2)
        subplot (2,2,1)
        plot(pi_c,FR)
        xlabel('pi_c')
        ylabel('FR')
        grid
        
        subplot (2,2,2)
        plot(pi_c,n0)
        xlabel('pi_c')
        ylabel('n0')
        grid
        
        subplot(2,2,3)
        plot(pi_c,nP)
        xlabel('pi_c')
        ylabel('nP')
        grid
        
    end  
end
%----------------------------------------------------------------------%

