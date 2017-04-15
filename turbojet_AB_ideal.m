% Turbojet ideal com afterburner
% aqui corremos em pi_c

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
Tt4  = 1600;            %temperaturas no fim da camara de combustao |K|
Tt7  = 2200;            %temperatura no fim do afterburner |K|
m0   = 0:0.5:3;         %M0 CONJUNTO FECHADO
M0   = m0';             %transposta
pi_c = 2:1:40;          %varredura em pi_c

j = length(M0);%cada linha do gr?fico CONJUNTO FECHADO
n = length(pi_c);
%-------------------------------------------------------------------------%

%saida--------------------------------------------------------------------%
%       empuxo = |N/(kg/s)|
%       f      = |adimensional| razao combustivel/ air
%       s      = |(mg/s)/N| consumo espef?cico de combust?vel por empuxo
%       nT     = |adimensional| eficiencia termica
%       nP     = |adimensional| eficiencia propulsiva
%       n0     = |adimensional| eficiencia total
%-------------------------------------------------------------------------%

%pre-alocando matrizes com zeros para preenche-las------------------------%
tal_r   = zeros(j,n); 
tal_lambidaAB = zeros(j,n); 
raz_1AB   = zeros(j,n);
raz_2AB   = zeros(j,n);
ftotal  = zeros(j,n);
sAB     = zeros(j,n);
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
    tal_lambidaAB = Tt7/T0;
    %raz_1 = V9/a0 AB
    raz_1AB     = sqrt((2/(y-1)).*tal_lambidaAB.*(1-(((tal_lambida)./(tal_r.*tal_c))./(tal_lambida-tal_r.*(tal_c-1)))));
    %raz_2 = F/m0ponto EMPUXO
    raz_2AB     = a0.*(raz_1AB-M0);
    f           = ((cp*T0)/hpr).*(tal_lambida-tal_r.*tal_c)/1000;
    ftotal      = ((cp*T0)/hpr).*(tal_lambidaAB-tal_r)/1000;
    fAB         = ftotal-f;
    sAB         = (ftotal./(raz_2AB)).*1000000;
    nT          = ((y-1).*cp.*T0.*(raz_1AB.^2-M0.^2))./(2.*ftotal.*hpr);
    nP          = (2.*M0)./(raz_1AB+M0);
    n0          = (nT.*nP);
    
    %Graficos---------------------------------------------------------%
        figure(1)
        subplot(2,2,1)
        plot(pi_c,raz_2AB)
        xlabel('pi_c')
        ylabel('empuxo')
        legend('M0_1','M0_2','M0_3')
        
        subplot(2,2,2)
        plot(pi_c,fAB,'r',pi_c,f,'b')
        xlabel('pi_c')
        ylabel('f')
       
        subplot(2,2,3)
        plot(pi_c,sAB)
        xlabel('pi_c')
        ylabel('s')
        
        subplot(2,2,4)
        plot(pi_c,nT,pi_c,nP,pi_c,n0)
        xlabel('pi_c')
        ylabel('n')
        
    end  
end
%----------------------------------------------------------------------%

