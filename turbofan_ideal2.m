% Turbofan ideal 
% aqui corremos em pi_f

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
Tt4  = 1670;            %temperaturas no fim da camara de combustao |K|
pi_f = 1:0.1:6;                
M0   = 0.9;             %M0 
pi_c = 24;          %varredura em pi_c
alfa = [3;4;5;10;12];               %CONJUNTO FECHADO razao de bypass (alfalinha - razao otima)

j = length(alfa);%cada linha do gr?fico CONJUNTO FECHADO
n = length(pi_f);   
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
    
    alfaestrela = (1./(tal_r.*(tal_f-1))).*(tal_lambida-tal_r.*(tal_c-1)-((tal_lambida)/(tal_r.*tal_c))-0.25.*((tal_r.*tal_f-1).^0.5+(tal_r-1).^0.5).^2);
    
    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_f,raz_2)
        xlabel('pi_f')
        ylabel('empuxo')
        grid
        
        %subplot(2,2,2)
        %plot(pi_f,f)
        %xlabel('pi_f')
        %ylabel('f')
       
        subplot(2,2,2)
        plot(pi_f,s)
        xlabel('pi_f')
        ylabel('s')
        grid
        
        subplot(2,2,3)
        plot(pi_f,nP,pi_f,n0)
        xlabel('pi_f')
        ylabel('nP')
        grid
        
        subplot(2,2,4)%ele corta y de 0 a 5
        plot(pi_f,FR)
        xlabel('pi_f')
        ylabel('FR')
        grid
        
        figure(2)
        plot(pi_f,alfaestrela)
        xlabel('pi_f')
        ylabel('Bypass otimo')
        grid
     
    end  
end
%----------------------------------------------------------------------%

