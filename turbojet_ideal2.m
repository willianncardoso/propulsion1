% Turbojet Ideal
% aqui corremos em pi_c

%zerando------------------------------------------------------------------%
clear variables
close all
clc
format long
%-------------------------------------------------------------------------%

% Dados iniciais conhecidos-----------------------------------------------%
T0   = 220;                   %temperatura inicial    |K|
y    = 1.4;                     %                       |J/kgK|
cp   = 1004;                    %                       |kJ/kg| 
hpr  = 42800;                              
Tt4  = 1400;       %temperaturas no fim da camara de combustao |K|
m0   = 0:0.5:3;    %M0 CONJUNTO FECHADO ()
M0   = m0';        %transposta  
pi_c = 0:0.5:20;   %varredura em pi_c - entra nas colunas das matrizes

j = length(M0);%cada linha do gr?fico CONJUNTO FECHADO
n = length(pi_c);
%-------------------------------------------------------------------------%

%saida--------------------------------------------------------------------%
%       empuxo = |N/(kg/s)|
%       f      = |adimensional| razao combustivel/ air
%       s      = |(mg/s)/N| consumo espeficico de combustivel por empuxo
%       nT     = |adimensional| eficiencia termica
%       nP     = |adimensional| eficiencia propulsiva
%       n0     = |adimensional| eficiencia total
%-------------------------------------------------------------------------%

%pre-alocando matrizes com zeros para preenche-las------------------------%
tal_r   = zeros(j,n); 
tal_lambida = zeros(j,n); 
raz_1   = zeros(j,n);
raz_2   = zeros(j,n);
f       = zeros(j,n);
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
    tal_r       = 1+((y-1)/2).*((M0).^2);
    tal_lambida = (Tt4)/(T0);
    tal_c       = (pi_c).^((y-1)/y);
    tal_t       = 1-(tal_r./tal_lambida).*(tal_c-1);
    %raz_1 = V9/a0 
    raz_1       = sqrt((2/(y-1)).*(tal_lambida./(tal_r*tal_c).*(tal_r.*tal_c.*tal_t-1)));
    %raz_2 = F/m0ponto EMPUXO
    raz_2       = a0.*(raz_1-M0);
    f           = ((cp*T0)/hpr).*(tal_lambida-tal_r.*tal_c)/1000;
    s           = (f./(raz_2)).*1000000;
    nT          = 1-(1./(tal_r.*tal_c));
    nP          = (2.*M0)./(raz_1+M0);
    n0          = (nT.*nP);
    
    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_c,raz_2)
        xlabel('pi_c')
        ylabel('empuxo')
        grid
        legend('0','0.5','1','1.5','2','2.5','3')
        
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
        plot(pi_c,nT,pi_c,nP,pi_c,n0)
        xlabel('pi_c')
        ylabel('n')
        grid
        legend('nT','nP','n0')
    
    end 
end
%----------------------------------------------------------------------%

